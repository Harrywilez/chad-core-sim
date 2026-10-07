"""The grading function: design → simulation → metrics → score.

Everything in the score comes from simulation values of *this* design; nothing is
fitted to the catalogue.  One evaluation (reactor preset, 1 m ball; the current is set so
that B_rms over the region of interest is `B_target` = 2 T, so every design is judged at
the same field strength and the current that takes is what buildability judges):

  1. exact Biot–Savart grid (25³, numba) and conductor distance map
  2. field-line topology on 96 seeds: closed fraction, mirror ratio, adiabatic prediction
  3. full-orbit Boris run of 15 keV D⁺ markers seeded uniformly in the vessel
     (r ∈ [0.18, 0.58], ≥ 0.03 from any conductor), isotropic pitch; 40 transits, extended
     to 80 and 120 while ≥ 10 % is still confined and the tail is a large share of τ
  4. two confinement times, in wall transits:
        τ_all  from all markers,
        τ_pass from the *passing* markers — those inside the loss cone of their own field
               line (ξ0² ≥ 1 − 1/R_line) plus every marker on a closed line: the population
               mirror trapping cannot hold.  A plasma needs them held too (they are refilled
               by collisions), so a trap that keeps only trapped pitch angles (a mirror, a
               dipole-like fat torus) is not rewarded as if it held everything.
     Each is  ∫₀ᵀ S dt + tail,  tail = S(T)/λ_tail  with λ_tail the loss rate over the last
     third of the run (floored at one expected loss) and the tail capped at S(T)·T, so an
     extrapolation can never exceed one more run-length.
        τ_c = √(τ_all · τ_pass)
  5. buildability = min(1, J_REBCO / J),  J = I_max / (π a²),  a = min(design wire radius,
     0.45 × clearance) — the same conductor radius the field model uses — and I_max includes
     the largest relative component current.  Invalid (score 0): conductors touching
     (clearance < 2 × wire radius) or outside the unit ball.

score = τ_c × buildability.  Multiply transits by R_wall / v for seconds (15 keV D⁺ at 1 m:
0.68 µs per transit).  Optional (toroidal designs): `surfaces=True` adds the Poincaré /
rotational-transform analysis of rung 1 as extra metrics.

Known limits (see agent/PROTOCOL.md): collisionless, vacuum field, no plasma pressure, no
self-consistent E; mirrors are judged harshly because their passing particles are lost by
construction — that is the physics, not a bug, but a collisional mirror model would be the
fair extension.
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Optional

import numpy as np

from .constants import DEUTERON, E_CHARGE
from .design import build, design_hash, validate
from .fieldlines import analyse, summarise
from .fields import FieldGrid
from .forces import CONDUCTOR_LIMITS_A_M2
from .particles import monoenergetic_ensemble, run_orbits, sample_positions_in_vessel
from .presets import PRESETS
from .runner import make_field, roi_mask

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "cache" / "designs"
LEDGER = ROOT / "results" / "design_ledger.jsonl"


def confinement_time(loss_times: np.ndarray, T: float, n: int) -> dict:
    """τ_c in the units of loss_times / T (censored integral + fitted exponential tail)."""
    lt = np.asarray(loss_times, dtype=float)
    finite = np.isfinite(lt)
    integral = float(np.sum(np.minimum(np.where(finite, lt, T), T))) / n          # ∫ S dt
    S_T = float(np.mean(~finite | (lt > T)))
    t0 = 2.0 * T / 3.0
    alive_t0 = int(np.sum(~finite | (lt > t0)))
    lost_tail = int(np.sum(finite & (lt > t0) & (lt <= T)))
    # exponential tail rate over the last third; floor at one expected loss among those alive
    lam = max(lost_tail, 1.0) / max(alive_t0, 1) / (T - t0) if alive_t0 > 0 else float("inf")
    tail = min(S_T / lam, S_T * T) if S_T > 0 else 0.0          # never more than one more run-length
    return {"tau_c": integral + tail, "integral": integral, "tail": tail, "S_T": S_T, "lambda_tail": lam,
            "tail_is_bound": lost_tail == 0, "tail_capped": bool(S_T > 0 and S_T / lam > S_T * T)}


def topology_batched(g: FieldGrid, seeds: np.ndarray, wall: float = 0.82, step: float = 0.008, max_length: float = 30.0, wire_clear: float = 0.008) -> dict:
    """Batched (all seeds at once) version of fieldlines.analyse on the grid: end kinds both
    ways, mirror ratio min(B_max+, B_max−)/B_seed, adiabatic retention √(1 − 1/R)."""
    n = len(seeds)
    B0 = g.Bmag_at(seeds)
    bmax = np.zeros((2, n)); ends = np.zeros((2, n), dtype=np.int8)   # 0 length, 1 wall, 2 wire, 3 grid-exit
    lengths = np.zeros((2, n))
    for di, direction in enumerate((1.0, -1.0)):
        x = seeds.copy(); alive = np.ones(n, bool); bm = B0.copy(); L = np.zeros(n)
        for _ in range(int(max_length / step)):
            idx = np.flatnonzero(alive)
            if not idx.size:
                break
            p = x[idx]

            def tang(q):
                b = g.B_at(q); nb = np.linalg.norm(b, axis=1); return direction * b / np.maximum(nb, 1e-300)[:, None]
            k1 = tang(p); k2 = tang(p + 0.5 * step * k1); k3 = tang(p + 0.5 * step * k2); k4 = tang(p + step * k3)
            pn = p + (step / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
            L[idx] += step
            bm[idx] = np.maximum(bm[idx], g.Bmag_at(pn))
            r = np.linalg.norm(pn, axis=1)
            hit_wall = r >= wall
            out = ~np.all(np.abs(pn) < g.half, axis=1)
            hit_wire = g.wire_distance_at(np.clip(pn, -g.half * 0.999, g.half * 0.999)) <= wire_clear
            x[idx] = pn
            for kind, m in ((1, hit_wall), (3, out), (2, hit_wire)):
                sel = idx[m & (ends[di, idx] == 0)]
                ends[di, sel] = kind; alive[sel] = False
        bmax[di] = bm; lengths[di] = L
    both_wall = (ends[0] == 1) & (ends[1] == 1)
    closed_like = (ends[0] == 0) & (ends[1] == 0)
    R = np.where(closed_like, np.inf, np.minimum(bmax[0], bmax[1]) / np.maximum(B0, 1e-300))
    pred = np.where(closed_like, 1.0, np.where(R > 1, np.sqrt(np.maximum(1 - 1 / np.maximum(R, 1e-300), 0)), 0.0))
    finite = R[np.isfinite(R)]
    return {"n_seeds": n, "frac_wall_both": float(both_wall.mean()), "frac_wire_any": float(((ends[0] == 2) | (ends[1] == 2)).mean()),
            "frac_closed": float(closed_like.mean()), "connection_length_median": float(np.median(lengths.sum(axis=0))),
            "mirror_ratio_median": float(np.median(finite)) if finite.size else float("nan"),
            "mirror_ratio_p90": float(np.percentile(finite, 90)) if finite.size else float("nan"),
            "predicted_adiabatic_retention": float(pred.mean()), "_R": R, "_closed": closed_like}


def evaluate(design: dict, transits: int = 20, n_particles: int = 160, grid_n: int = 25, energy_eV: float = 15e3,
             surfaces: bool = False, seed: int = 12345, verbose: bool = True, log: bool = True, B_target: float = 2.0,
             max_transits: int = 120, position_seed: int = 271828, topology_seed: int = 31415) -> dict:
    """Evaluate with independent velocity, particle-position, and topology seeds.

    Defaults preserve the historical initial conditions. Change all three seeds
    for held-out validation; changing ``seed`` alone changes velocities only.
    """
    t_start = time.time()
    v = validate(design)
    out = {"name": design.get("name", "design"), "hash": v["hash"], "design": design, "validation": v, "transits": transits,
           "n_particles": n_particles, "grid_n": grid_n, "energy_eV": energy_eV, "B_target_T": B_target,
           "seed": seed, "position_seed": position_seed, "topology_seed": topology_seed,
           "max_transits": max_transits, "surfaces_requested": surfaces, "evaluation_version": "2.3"}
    if not v["valid"]:
        out.update({"score": 0.0, "tau_c": 0.0, "invalid": True, "seconds": time.time() - t_start})
        if log:
            _log(out)
        return out
    w = build(design)
    pr = PRESETS["reactor"]()
    CACHE.mkdir(parents=True, exist_ok=True)
    gpath = CACHE / f"grid_{v['hash']}_n{grid_n}.npz"
    if gpath.exists():
        d = np.load(gpath)
        g = FieldGrid.__new__(FieldGrid)
        g.winding, g.half, g.n = w, 0.87, grid_n
        g.axis = np.linspace(-0.87, 0.87, grid_n); g.dx = float(g.axis[1] - g.axis[0])
        xx, yy, zz = np.meshgrid(g.axis, g.axis, g.axis, indexing="ij")
        g.points = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
        g.B = d["B"]; g.A = None; g.Bmag = np.linalg.norm(g.B, axis=3); g.wire_distance = d["wire_distance"]
    else:
        g = FieldGrid(w, 0.87, grid_n, with_A=False)
        np.savez_compressed(gpath, B=g.B, wire_distance=g.wire_distance)
    t_grid = time.time() - t_start
    # --- topology ---
    seed_pool = sample_positions_in_vessel(600, 0.18, 0.58, np.random.default_rng(topology_seed))
    seeds = seed_pool[g.wire_distance_at(seed_pool) > 0.03][:96]
    topo = topology_batched(g, seeds)
    # --- field statistics at reactor scale ---
    mask = roi_mask(g)
    Bm = g.Bmag.ravel()[mask]
    # Design radii are in unit-ball coordinates; collision and buildability
    # must use the same declared conductor, scaled into reactor metres.
    requested_radius_m = float(design.get("wire_radius", 0.005)) * pr.unit_ball_m
    field0, a = make_field(g, pr, w, pr.current_A, requested_radius_m)
    brms0 = float(abs(field0.gain) * np.sqrt(np.mean(Bm**2)))
    current = pr.current_A * B_target / max(brms0, 1e-12)          # current (per unit component current) for B_rms = B_target
    field, a = make_field(g, pr, w, current, requested_radius_m)
    brms = float(abs(field.gain) * np.sqrt(np.mean(Bm**2)))
    # --- particles ---
    pool = sample_positions_in_vessel(2000, 0.18, 0.58, np.random.default_rng(position_seed))
    pos = pool[g.wire_distance_at(pool) > 0.03][:n_particles]
    ens = monoenergetic_ensemble(DEUTERON, n_particles, energy_eV, pos, np.random.default_rng(seed))
    vD = math.sqrt(2 * energy_eV * E_CHARGE / DEUTERON.mass_kg)
    transit = 0.82 / vD
    omega_ref = E_CHARGE / DEUTERON.mass_kg * max(brms, 1e-12)
    dt = min(pr.dt_s, 0.15 / omega_ref)
    t_max = transits * transit
    try:
        from .fastorbits import HAVE_NUMBA, fast_orbits
    except Exception:  # pragma: no cover
        HAVE_NUMBA = False
    out["orbit_backend"] = "numba" if HAVE_NUMBA else "numpy"
    # adaptive length: extend the run (20 → 40 → 60 … up to max_transits) while the extrapolated tail
    # exceeds what was actually integrated, so a flat stretch in a short run cannot inflate τ_c
    T_run = max(transits, 40)
    while True:
        ens_run = monoenergetic_ensemble(DEUTERON, n_particles, energy_eV, pos, np.random.default_rng(seed))
        if HAVE_NUMBA:
            lt_s, kinds = fast_orbits(ens_run, field, T_run * transit, dt, 0.82)
            lt = lt_s / transit
            counts = {k: int(np.sum(kinds == k)) for k in ("retained", "wall", "wire", "grid-exit", "numerical")}
            alive = ~np.isfinite(lt_s)
            energy_drift = float("nan")
        else:
            res = run_orbits(ens_run, field, T_run * transit, dt, 0.82)
            # Keep original marker IDs: the passing mask below is marker-indexed.
            # res.loss_times contains only lost particles and cannot be padded safely.
            lt = ens_run.loss_time.copy() / transit
            counts, alive, energy_drift = res.counts, ens_run.alive, res.energy_drift_rel_max
        tau = confinement_time(lt, float(T_run), n_particles)
        # extend while a plateau is still unresolved: something is still confined and the tail is a large share
        if T_run >= max_transits or not alive.any() or (tau["S_T"] < 0.1 and tau["tail"] <= tau["integral"]):
            break
        if tau["S_T"] < 0.1 or tau["tail"] <= 0.5 * tau["integral"]:
            break
        T_run = min(T_run * 2, max_transits)
    ens = ens_run
    S = {f"S{c}": float(np.mean(~np.isfinite(lt) | (lt > c))) for c in (4, 8, 20, 40, 60, 80, 120) if c <= T_run}
    out["transits_run"] = T_run
    # "passing" markers: the ones mirror trapping cannot hold — inside the loss cone of their own field
    # line (ξ0² ≥ 1 − 1/R_line, R_line = min(B_max+, B_max−)/B0 from the topology tracer), and every
    # marker on a closed line (no loss cone there).  τ_pass is what closed lines / surfaces buy.
    ens0 = monoenergetic_ensemble(DEUTERON, n_particles, energy_eV, pos, np.random.default_rng(seed))
    B0 = field.B_at(ens0.x)
    b0 = B0 / np.maximum(np.linalg.norm(B0, axis=1), 1e-300)[:, None]
    xi0 = np.abs(np.einsum("ni,ni->n", ens0.v, b0)) / np.linalg.norm(ens0.v, axis=1)
    topo_m = topology_batched(g, pos)
    R_line = topo_m["_R"]
    passing = topo_m["_closed"] | (xi0**2 >= 1.0 - 1.0 / np.maximum(R_line, 1.0 + 1e-9))
    tau_pass = confinement_time(lt[passing], float(T_run), int(passing.sum())) if passing.sum() >= 8 else {"tau_c": tau["tau_c"], "note": "too few passing markers"}
    S_pass = {f"S{c}": float(np.mean(~np.isfinite(lt[passing]) | (lt[passing] > c))) for c in (4, 8, 20, 40, 60, 80, 120) if c <= T_run}
    core = [float(np.mean(np.linalg.norm(ens.x[alive], axis=1) < 0.35))] if alive.any() else []
    # --- buildability ---
    # `current` is the physical current of the FIRST winding (ScaledGridField normalises the grid by its
    # relative amplitude), so the physical current per unit relative amplitude is current / |rel_0| and the
    # largest physical current is that times the largest relative amplitude.  (Before v2.3 Imax was
    # current × rel_max, which double-counted a first component whose relative current was not 1 — a uniform
    # rescaling of all relative currents, physically a no-op, changed the score; events 122/126/139.)
    ref_rel = abs(float(w.windings[0].current_A))
    I_unit = current / ref_rel
    scaled = w.scaled(1.0).with_current(current).with_wire_radius(a)       # physical per-winding currents (relative amplitudes kept)
    rel_max = max(abs(float(wd.current_A)) for wd in w.windings)
    Imax = I_unit * rel_max
    J = Imax / (math.pi * a * a)
    # ROI diagnostics (not in the score): the 2 T normalisation averages |B| over grid points ≥ 0.08 from any
    # conductor, so an auxiliary coil parked in a weak-field region removes weak points from the average and
    # lowers the current needed for B_target — record how much of the shell the mask keeps and a mask-free median.
    r_all = np.linalg.norm(g.points, axis=1)
    shell = (r_all >= 0.20) & (r_all <= 0.65)
    B_median_shell = float(abs(field.gain) * np.median(g.Bmag.ravel()[shell]))
    roi_fraction = float(mask.sum() / max(1, shell.sum()))
    J_lim = CONDUCTOR_LIMITS_A_M2["REBCO tape at 20 T, 20 K (engineering)"]
    build_factor = min(1.0, J_lim / J)
    tau_c = math.sqrt(max(tau["tau_c"], 0.0) * max(tau_pass["tau_c"], 0.0))
    score = tau_c * build_factor
    out.update({
        "topology": {k: val for k, val in topo.items() if not k.startswith("_")},
        "field": {"B_rms_roi_T": brms, "B_rms_per_364kA_T": brms0, "current_A": I_unit, "current_first_winding_A": current,
                  "reference_relative_current": ref_rel, "B_max_over_rms": float(Bm.max() / np.sqrt(np.mean(Bm**2))),
                  "B_median_shell_T": B_median_shell, "roi_fraction": roi_fraction, "roi_points": int(mask.sum()),
                  "wire_radius_used_m": a, "J_A_m2": J, "joule_power_W_copper": scaled.joule_power_W()},
        "particles": {"survival": S, "survival_passing": S_pass, "n_passing": int(passing.sum()), "tau_all": tau["tau_c"], "tau_passing": tau_pass["tau_c"],
                      "tau_c_transits": tau_c, "tau_detail": tau, "tau_passing_detail": tau_pass, "core_fraction_final": float(np.mean(core)) if core else float("nan"),
                      "energy_drift_rel_max": energy_drift, "counts": counts, "transit_seconds_reactor": transit,
                      "loss_times_transits": [float(x) if np.isfinite(x) else None for x in lt]},
        "buildability": {"factor": build_factor, "clearance": v["clearance"], "wire_radius_used": a, "J_over_REBCO": J / J_lim, "I_max_A": Imax},
        "tau_c": tau_c, "score": score, "invalid": False, "grid_seconds": t_grid,
    })
    if surfaces:
        try:
            from .surfaces import analyse_from_axis, torus_frame
            R0, rt = torus_frame(w.active_points)
            sres = analyse_from_axis(w, R0, rt, n_seeds=16, max_turns=30)
            out["surfaces"] = {k: sres[k] for k in ("axis_R", "axis_z", "frac_surface", "frac_island_chaotic", "frac_open", "outermost_surface_r", "iota_axis", "iota_edge", "well_depth", "iota_profile")}
            # plasma-aware score (rung 3): closed-line configurations without enough rotational transform lose the
            # plasma by polarization on τ_pol; a design with surfaces and ι above the shorting threshold is exempt
            from .polarization import polarization_loss
            R_m, a_m = max(float(sres["axis_R"]), 0.05), max(float(sres["outermost_surface_r"]), 0.02)
            pol = polarization_loss(energy_eV, brms, R_m, a_m, DEUTERON.mass_kg, abs(float(sres["iota_edge"])) if sres["frac_surface"] > 0 else 0.0)
            has_surfaces = sres["frac_surface"] >= 0.2 and pol["shorted"]
            tau_pol_wall_transits = pol["tau_pol_s"] / transit
            # v2.2: the cap applies to every toroidal-field design without surfaces — closed lines (frac_closed ≥ 0.2)
            # OR a leaky toroidal solenoid whose lines drift out slowly (mirror ratio ≈ 1: nothing is mirror-trapped, the
            # passing markers are only held by the slow leak, and a plasma polarizes on τ_pol regardless).  Mirror-type
            # designs (mirror ratio ≥ 1.5) get no factor: their passing-particle loss is already in τ_pass.
            toroidal_like = topo["frac_closed"] >= 0.2 or topo["mirror_ratio_median"] < 1.5
            factor = 1.0 if has_surfaces else (min(1.0, tau_pol_wall_transits / max(tau_c, 1e-9)) if toroidal_like else 1.0)
            out["plasma"] = {"polarization": pol, "tau_pol_wall_transits": tau_pol_wall_transits, "surfaces_ok": bool(has_surfaces), "factor": factor,
                             "toroidal_like": bool(toroidal_like)}
            out["score_plasma"] = score * factor
        except Exception as e:  # pragma: no cover
            out["surfaces"] = {"error": str(e)}
    out["seconds"] = time.time() - t_start
    if verbose:
        print(f"{out['name'][:44]:44s} τ_c={tau_c:6.1f} (all {tau['tau_c']:.0f} · pass {tau_pass['tau_c']:.0f}; {T_run}T) score={score:6.1f}  " + " ".join(f"{k}={val:.2f}" for k, val in S.items() if k in ('S4', 'S20', 'S120'))
              + f"  closed={topo['frac_closed']:.2f} pred={topo['predicted_adiabatic_retention']:.2f} I={I_unit/1e3:.0f}kA build={build_factor:.2f} ({out['seconds']:.0f}s)", flush=True)
    if log:
        _log(out)
    return out


def _log(rec: dict):
    LEDGER.parent.mkdir(parents=True, exist_ok=True)
    slim = {k: v for k, v in rec.items() if k not in ("particles",)}
    if "particles" in rec:
        slim["particles"] = {k: v for k, v in rec["particles"].items() if k != "loss_times_transits"}
    slim["time"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
    with open(LEDGER, "a") as fh:
        fh.write(json.dumps(slim, default=float) + "\n")


def main():  # CLI: python -m ccsim.evaluate design.json [--transits 20] [--particles 160] [--surfaces]
    import argparse
    ap = argparse.ArgumentParser(description="Evaluate a design JSON and print its metrics/score.")
    ap.add_argument("design", help="path to a design JSON file (or '-' for stdin)")
    ap.add_argument("--transits", type=int, default=20)
    ap.add_argument("--max-transits", type=int, default=120)
    ap.add_argument("--particles", type=int, default=160)
    ap.add_argument("--grid", type=int, default=25)
    ap.add_argument("--surfaces", action="store_true")
    ap.add_argument("--seed", type=int, default=12345, help="velocity seed")
    ap.add_argument("--position-seed", type=int, default=271828, help="particle-position seed")
    ap.add_argument("--topology-seed", type=int, default=31415, help="field-line seed positions")
    ap.add_argument("--out", help="write the full result JSON here")
    a = ap.parse_args()
    import sys
    design = json.load(sys.stdin if a.design == "-" else open(a.design))
    res = evaluate(design, transits=a.transits, n_particles=a.particles, grid_n=a.grid, surfaces=a.surfaces,
                   max_transits=a.max_transits, seed=a.seed, position_seed=a.position_seed, topology_seed=a.topology_seed)
    if a.out:
        json.dump(res, open(a.out, "w"), indent=1, default=float)
    summary = {k: res[k] for k in ("name", "hash", "score", "tau_c", "invalid")} | {"survival": res.get("particles", {}).get("survival"), "validation": res["validation"]}
    if "score_plasma" in res:
        summary["score_plasma"] = res["score_plasma"]; summary["surfaces"] = {k: res["surfaces"].get(k) for k in ("frac_surface", "outermost_surface_r", "iota_axis", "iota_edge", "well_depth")}
    print(json.dumps(summary, default=float))


if __name__ == "__main__":
    main()
