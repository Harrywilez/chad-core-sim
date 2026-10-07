"""exp09 — follow-ups on the exp08 playbook.

  a) energy scan (5 / 15 / 45 keV D⁺) for the three best families: drift-limited
     confinement should scale like 1/ρ_L (time to drift across the trap), loss-cone
     limited confinement should not care.
  b) coarse-mesh mirror pairs (4 / 6 / 8 circuits, opposed and unequal currents) —
     the woven pair with few circuits is a two-strand helical torus; does the
     discrete helical ripple give closed flux surfaces with rotational transform?
  c) rotational transform ι = Δθ_pol / Δφ_tor measured on long field lines for the
     toroidal configurations (torus centre circle from the strand's ρ range).
  d) 60-transit survival for the best configuration.

Usage: python3 experiments/exp09_followups.py [a|b|all]
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ccsim.constants import DEUTERON, E_CHARGE
from ccsim.fieldlines import analyse, summarise
from ccsim.fields import FieldGrid
from ccsim.geometry import (MultiWinding, Winding, assembly_windings, close_with_return, codex_toroidal, hopf_continued,
                            hopf_mirror_pair, hopf_mirror_windings, hopf_torus)
from ccsim.particles import monoenergetic_ensemble, run_orbits, sample_positions_in_vessel
from ccsim.presets import PRESETS
from ccsim.runner import make_field, roi_mask

R = ROOT / "results"
C = ROOT / "cache"
which = sys.argv[1] if len(sys.argv) > 1 else "all"
GRID_N = 25
pr = PRESETS["reactor"]()
pool = sample_positions_in_vessel(2000, 0.18, 0.58, np.random.default_rng(271828))
seed_pool = sample_positions_in_vessel(600, 0.18, 0.58, np.random.default_rng(31415))


def single(name, pts, meta):
    closed, n_active = close_with_return(pts, 3.8)
    return Winding(name, closed, 1.0, 0.005, closed=True, active_count=n_active, meta=meta)


def mirror_unequal(name, alpha, **kw):
    """Original strand at +I, mirror strand at alpha·I."""
    pair = hopf_mirror_pair(**kw)
    ws = []
    for k, (pts, sgn) in enumerate(zip(pair, (1.0, alpha))):
        closed, n_active = close_with_return(pts, 3.8)
        ws.append(Winding(f"{name} #{k}", closed, sgn, 0.005, closed=True, active_count=n_active, meta={"family": "hopf-mirror", "index": k, "sign": sgn}))
    return MultiWinding(name, ws, {"family": "hopf-mirror", "alpha": alpha, **kw})


def load_grid(w):
    key = "".join(c if c.isalnum() else "_" for c in w.name)[:70]
    path = C / f"grid_{key}_n{GRID_N}.npz"
    if path.exists():
        d = np.load(path)
        g = FieldGrid.__new__(FieldGrid)
        g.winding, g.half, g.n = w, 0.87, GRID_N
        g.axis = np.linspace(-0.87, 0.87, GRID_N); g.dx = float(g.axis[1] - g.axis[0])
        xx, yy, zz = np.meshgrid(g.axis, g.axis, g.axis, indexing="ij")
        g.points = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
        g.B = d["B"]; g.A = None; g.Bmag = np.linalg.norm(g.B, axis=3); g.wire_distance = d["wire_distance"]
        return g
    g = FieldGrid(w, 0.87, GRID_N, with_A=False)
    np.savez_compressed(path, B=g.B, wire_distance=g.wire_distance)
    return g


def particle_run(w, g, energy_eV, transits, n_part=160):
    field, a = make_field(g, pr, w, pr.current_A, pr.wire_radius_m)
    mask = roi_mask(g)
    brms = float(abs(field.gain) * np.sqrt(np.mean(g.Bmag.ravel()[mask] ** 2)))
    pos = pool[g.wire_distance_at(pool) > 0.03][:n_part]
    ens = monoenergetic_ensemble(DEUTERON, n_part, energy_eV, pos, np.random.default_rng(12345))
    v = math.sqrt(2 * energy_eV * E_CHARGE / DEUTERON.mass_kg)
    omega_ref = E_CHARGE / DEUTERON.mass_kg * max(brms, 1e-12)
    dt = min(pr.dt_s, 0.15 / omega_ref)
    t_max = transits * 0.82 / v
    res = run_orbits(ens, field, t_max, dt, 0.82)
    lost = res.loss_times
    S = {f"S{c}": float(1 - np.sum(lost <= c * 0.82 / v) / n_part) for c in (4, 8, 20, 40, 60) if c <= transits}
    med = float(np.median(lost[np.isfinite(lost)])) * v / 0.82 if np.isfinite(lost).any() else float("inf")
    rho = DEUTERON.mass_kg * v / (E_CHARGE * brms)
    return {"energy_eV": energy_eV, "transits": transits, "survival": S, "median_loss_transits": med, "rho_L_m": rho, "B_rms_T": brms,
            "mean_confined_transits": float(np.mean(np.minimum(np.where(np.isfinite(lost), lost, t_max), t_max)) * v / 0.82)}


def rotational_transform(g, w, n_lines=24):
    """ι on long field lines.  Torus centre circle radius from the strand's ρ-range (round torus)."""
    act = w.active_points
    rho = np.hypot(act[:, 0], act[:, 1])
    R0 = 0.5 * (rho.min() + rho.max())
    seeds = seed_pool[g.wire_distance_at(seed_pool) > 0.03][:n_lines]
    out = []
    for s in seeds:
        pts = [np.asarray(s, dtype=float)]
        p = pts[0].copy(); ended = "length"
        for _ in range(6000):  # 36 ball units of line
            def tangent(q):
                b = g.B_at(q[None, :])[0]; n = np.linalg.norm(b)
                return b / n if n > 0 else np.zeros(3)
            k1 = tangent(p); k2 = tangent(p + 0.003 * k1); k3 = tangent(p + 0.003 * k2); k4 = tangent(p + 0.006 * k3)
            p = p + 0.001 * (k1 + 2 * k2 + 2 * k3 + k4)
            if np.linalg.norm(p) >= 0.82: ended = "wall"; break
            if not np.all(np.abs(p) < g.half): ended = "grid-exit"; break
            if g.wire_distance_at(p[None, :])[0] <= 0.008: ended = "wire"; break
            pts.append(p.copy())
        pts = np.array(pts)
        if len(pts) < 200:
            continue
        phi = np.unwrap(np.arctan2(pts[:, 1], pts[:, 0]))
        th = np.unwrap(np.arctan2(pts[:, 2], np.hypot(pts[:, 0], pts[:, 1]) - R0))
        if abs(phi[-1] - phi[0]) < 2 * np.pi:
            continue
        out.append({"iota": float((th[-1] - th[0]) / (phi[-1] - phi[0])), "toroidal_turns": float((phi[-1] - phi[0]) / (2 * np.pi)),
                    "ended": ended, "length": float(len(pts) * 0.006)})
    if not out:
        return {"n_lines": 0}
    io = np.array([o["iota"] for o in out])
    return {"n_lines": len(out), "iota_median": float(np.median(io)), "iota_p10": float(np.percentile(io, 10)), "iota_p90": float(np.percentile(io, 90)),
            "R0": float(R0), "turns_median": float(np.median([o["toroidal_turns"] for o in out]))}


out_path = R / f"exp09_{which}.json"
results = json.load(open(out_path)) if out_path.exists() else {}

jobs = []
if which in ("a", "all"):
    best = {
        "hopf mirror woven η=0.70 24×1 opposed": hopf_mirror_windings(1.0, 0.005, 1.0, -1, name="hopf mirror woven η=0.70 24×1 opposed", eta=0.70, circuits=24, mode="woven"),
        "precess × pair-side-crossed, core 0.45": MultiWinding("precess × pair-side-crossed, core 0.45",
                                                              assembly_windings(codex_toroidal("precess", 0.22, 24.0, 3.0), "pair-side-crossed", 1.0, 0.005, 1.0, "precess"),
                                                              {"family": "assembly", "core": "precess", "layout": "pair-side-crossed"}),
        "hopf continued 720° / 48 circuits": single("hopf continued 720° / 48 circuits", hopf_continued(0.82, 0.64, 720, 48), {"family": "hopf-continued"}),
    }
    for name, w in best.items():
        for E in (5e3, 15e3, 45e3):
            jobs.append(("energy", name, w, E, 20))
    jobs.append(("long", "hopf mirror woven η=0.70 24×1 opposed", best["hopf mirror woven η=0.70 24×1 opposed"], 15e3, 60))
if which in ("b", "all"):
    for N in (4, 6, 8):
        jobs.append(("mesh", f"hopf mirror woven η=0.70 {N}×1 opposed", hopf_mirror_windings(1.0, 0.005, 1.0, -1, name=f"hopf mirror woven η=0.70 {N}×1 opposed", eta=0.70, circuits=N, mode="woven", delta_eta=0.05), 15e3, 20))
    for alpha in (-0.5, 0.0):
        jobs.append(("mesh", f"hopf mirror woven η=0.70 24×1, mirror at {alpha:+.1f} I", mirror_unequal(f"hopf mirror woven η=0.70 24×1, mirror at {alpha:+.1f} I", alpha, eta=0.70, circuits=24, mode="woven"), 15e3, 20))
    jobs.append(("mesh", "hopf mirror woven η=0.70 6×1, mirror at -0.5 I", mirror_unequal("hopf mirror woven η=0.70 6×1, mirror at -0.5 I", -0.5, eta=0.70, circuits=6, mode="woven", delta_eta=0.05), 15e3, 20))

for kind, name, w, E, transits in jobs:
    key = f"{name} | {E/1e3:.0f} keV | {transits} transits"
    if key in results:
        continue
    t0 = time.time()
    g = load_grid(w)
    rec = {"kind": kind, "name": name, "clearance": w.clearance_m()}
    if kind == "mesh" or (kind == "energy" and E == 15e3):
        seeds = seed_pool[g.wire_distance_at(seed_pool) > 0.03][:96]
        topo = summarise(analyse(g, seeds, wall_radius=0.82, step=0.006, max_length=40.0, wire_clear=0.008))
        rec["topology"] = {k: v for k, v in topo.items() if not k.startswith("per_seed")}
        if "mirror" in name or "torus" in name:
            rec["transform"] = rotational_transform(g, w)
    rec["run"] = particle_run(w, g, E, transits)
    rec["seconds"] = time.time() - t0
    results[key] = rec
    S = rec["run"]["survival"]
    print(f"{key:70s} " + " ".join(f"{k}={v:.2f}" for k, v in S.items()) + f" med={rec['run']['median_loss_transits']:.1f} ρ={rec['run']['rho_L_m']*1e3:.1f} mm"
          + (f" closed={rec['topology']['frac_closed']:.2f}" if "topology" in rec else "") + (f" ι={rec['transform'].get('iota_median', float('nan')):.3f} (n={rec['transform']['n_lines']})" if "transform" in rec else "")
          + f" ({time.time()-t0:.0f}s)", flush=True)
    json.dump(results, open(out_path, "w"), default=float)
print("done", which)
