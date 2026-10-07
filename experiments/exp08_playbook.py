"""exp08 — 'play with it': parameter sweeps over the Hopf family, the precessing-loop
assemblies, the new Hopf mirror pair (meshed torus) and the new geometries.

Every configuration gets the same treatment (unit ball, wire radius 0.005, 25³ grid):
  * field-line topology on 96 seeds (closed / wall / mirror ratio / adiabatic prediction)
  * reactor-scale D⁺ (15 keV, I = 364 kA per circuit) full-orbit run, 160 particles,
    20 wall transits, survival S4 / S8 / S20, median loss time, energy drift,
    time-averaged core concentration (fraction of live particles inside r < 0.35)
  * B_rms in the ROI, B_max/B_rms (field uniformity), Joule power at reactor scale

Usage: python3 experiments/exp08_playbook.py <part>   (part = 0/1 splits the list; 'all' runs everything)
Results accumulate in results/exp08_<part>.json and are merged by exp08_collect.py.
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
from ccsim.geometry import (ASSEMBLY_LAYOUTS, MultiWinding, Winding, assembly_windings, baseball_seam, close_with_return,
                            codex_toroidal, hopf_continued, hopf_drift, hopf_link, hopf_mirror_windings, hopf_torus,
                            multi_windings, picket_fence, solenoid, sphere_winding, yin_yang)
from ccsim.particles import monoenergetic_ensemble, run_orbits, sample_positions_in_vessel
from ccsim.presets import PRESETS
from ccsim.runner import make_field, roi_mask

R = ROOT / "results"
C = ROOT / "cache"
part = sys.argv[1] if len(sys.argv) > 1 else "all"
GRID_N = 25


def single(name, pts, meta):
    closed, n_active = close_with_return(pts, 3.8)
    return Winding(name, closed, 1.0, 0.005, closed=True, active_count=n_active, meta=meta)


def multi(name, parts, family, meta=None):
    return multi_windings(parts, 1.0, 0.005, 1.0, name, family, meta=meta)


configs = []
# --- A. Hopf mirror pair (meshed torus) ---------------------------------------------------
for mode in ("woven", "nested"):
    for sense in (1, -1):
        configs.append(hopf_mirror_windings(1.0, 0.005, 1.0, sense, name=f"hopf mirror {mode} η=0.70 24×1 {'same' if sense > 0 else 'opposed'}",
                                            eta=0.70, circuits=24, revolutions=1.0, mode=mode))
for eta in (0.60, 0.80):
    configs.append(hopf_mirror_windings(1.0, 0.005, 1.0, -1, name=f"hopf mirror woven η={eta:.2f} 24×1 opposed", eta=eta, circuits=24, mode="woven"))
    configs.append(hopf_mirror_windings(1.0, 0.005, 1.0, 1, name=f"hopf mirror woven η={eta:.2f} 24×1 same", eta=eta, circuits=24, mode="woven"))
for circ in (16, 36):
    configs.append(hopf_mirror_windings(1.0, 0.005, 1.0, -1, name=f"hopf mirror woven η=0.70 {circ}×1 opposed", eta=0.70, circuits=circ, mode="woven"))
# --- B. Hopf torus η sweep ------------------------------------------------------------------
for eta in (0.55, 0.62, 0.70, 0.78):
    configs.append(single(f"hopf torus η={eta:.2f}, 24×1", hopf_torus(eta, 24, 1.0), {"family": "hopf-torus", "eta": eta, "circuits": 24, "revolutions": 1}))
configs.append(single("hopf torus η=0.70, 25×3", hopf_torus(0.70, 25, 3.0), {"family": "hopf-torus", "eta": 0.70, "circuits": 25, "revolutions": 3}))
# --- C. Hopf continued sweep ----------------------------------------------------------------
for sweep, circ in ((360, 24), (540, 36), (720, 48), (720, 32), (720, 64)):
    configs.append(single(f"hopf continued {sweep}° / {circ} circuits", hopf_continued(0.82, 0.64, sweep, circ),
                          {"family": "hopf-continued", "sweep": sweep, "circuits": circ, "eta1": 0.64}))
for eta1 in (0.72, 0.56):
    configs.append(single(f"hopf continued 720° / 48, η→{eta1:.2f}", hopf_continued(0.82, eta1, 720, 48),
                          {"family": "hopf-continued", "sweep": 720, "circuits": 48, "eta1": eta1}))
# --- D. Hopf drift circuits ------------------------------------------------------------------
for circ in (6, 10, 16):
    configs.append(single(f"hopf drift Shallow 180° / {circ} circuits", hopf_drift("s64_180", circ), {"family": "hopf-drift", "variant": "s64_180", "circuits": circ}))
# --- E. Precessing-loop pair variants --------------------------------------------------------
base = codex_toroidal("precess", 0.22, 24.0, 3.0)
for cs in (0.35, 0.45, 0.55):
    ws = assembly_windings(base, "pair-side-crossed", 1.0, 0.005, 1.0, "precess", core_scale=cs)
    configs.append(MultiWinding(f"precess × pair-side-crossed, core {cs:.2f}", ws, {"family": "assembly", "core": "precess", "layout": "pair-side-crossed", "core_scale": cs}))
for inward, rot in ((0.15, 24.0), (0.30, 24.0), (0.22, 12.0), (0.22, 36.0)):
    core = codex_toroidal("precess", inward, rot, 3.0)
    ws = assembly_windings(core, "pair-side-crossed", 1.0, 0.005, 1.0, f"precess i{inward} r{rot}")
    configs.append(MultiWinding(f"precess(inward {inward:.2f}, rot {rot:.0f}°) × pair-side-crossed", ws,
                                {"family": "assembly", "core": "precess", "layout": "pair-side-crossed", "inward": inward, "rotation_deg": rot}))
ws = assembly_windings(base, "pair-axial-same", 1.0, 0.005, 1.0, "precess", core_scale=0.55)
configs.append(MultiWinding("precess × pair-axial-same, core 0.55", ws, {"family": "assembly", "core": "precess", "layout": "pair-axial-same", "core_scale": 0.55}))
# --- F. New geometries -----------------------------------------------------------------------
configs.append(single("baseball seam (6 turns, 40°)", baseball_seam(), {"family": "baseball", "turns": 6, "amplitude_deg": 40}))
configs.append(single("baseball seam (6 turns, 60°)", baseball_seam(amplitude_deg=60.0), {"family": "baseball", "turns": 6, "amplitude_deg": 60}))
configs.append(multi("yin-yang pair (4+4 turns)", yin_yang(), "yin-yang", {"turns": 4}))
configs.append(multi("picket fence 5 rings", picket_fence(5), "picket", {"rings": 5}))
configs.append(multi("picket fence 7 rings", picket_fence(7), "picket", {"rings": 7}))
configs.append(single("sphere winding 24 turns (uniform B)", sphere_winding(), {"family": "sphere", "turns": 24}))
configs.append(multi("Hopf link (two linked rings)", hopf_link(), "link"))
configs.append(multi("Hopf link, opposed", [(hopf_link()[0], 1), (hopf_link()[1], -1)], "link", {"sense": -1}))
bb = baseball_seam()
ws = assembly_windings(bb, "pair-side-crossed", 1.0, 0.005, 1.0, "baseball")
configs.append(MultiWinding("baseball × pair-side-crossed", ws, {"family": "assembly", "core": "baseball", "layout": "pair-side-crossed"}))

if part == "user":
    # the configuration John picked out in the Lab: yin-yang pair (0.76 / 0.60, 40°, 4 turns each) in the cube-8 assembly, core scale 0.5
    from ccsim.geometry import assembly
    def yy_cube(layout, inner_sign, label):
        outer, inner = yin_yang()
        parts = [(pts, sgn) for pts, sgn in assembly(outer, layout, 0.5)] + [(pts, sgn * inner_sign) for pts, sgn in assembly(inner, layout, 0.5)]
        return multi(label, parts, "assembly", {"core": "yin-yang", "layout": layout, "core_scale": 0.5, "inner_sign": inner_sign})
    configs = [yy_cube("cube-8", 1, "yin-yang × cube-8 (John's pick)"), yy_cube("cube-8", -1, "yin-yang (inner opposed) × cube-8"),
               yy_cube("cube-8-checker", 1, "yin-yang × cube-8-checker"), yy_cube("pair-side-crossed", 1, "yin-yang × pair-side-crossed")]
elif part != "all":
    k = int(part)
    configs = [c for i, c in enumerate(configs) if i % 2 == k]

out_path = R / f"exp08_{part}.json"
results = json.load(open(out_path)) if out_path.exists() else {}
pool = sample_positions_in_vessel(2000, 0.18, 0.58, np.random.default_rng(271828))
seed_pool = sample_positions_in_vessel(600, 0.18, 0.58, np.random.default_rng(31415))
pr = PRESETS["reactor"]()
v_D = math.sqrt(2 * 15e3 * E_CHARGE / DEUTERON.mass_kg)
N_PART = 160
TRANSITS = 20

for w in configs:
    if w.name in results:
        continue
    t0 = time.time()
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
    else:
        g = FieldGrid(w, 0.87, GRID_N, with_A=False)
        np.savez_compressed(path, B=g.B, wire_distance=g.wire_distance)
    t_grid = time.time() - t0
    # topology
    seeds = seed_pool[g.wire_distance_at(seed_pool) > 0.03][:96]
    lines = analyse(g, seeds, wall_radius=0.82, step=0.006, max_length=40.0, wire_clear=0.008)
    topo = summarise(lines)
    # field statistics
    mask = roi_mask(g)
    Bm = g.Bmag.ravel()[mask]
    field, a = make_field(g, pr, w, pr.current_A, pr.wire_radius_m)
    brms = float(abs(field.gain) * np.sqrt(np.mean(Bm**2)))
    bstats = {"B_rms_roi_T": brms, "B_max_over_rms": float(Bm.max() / np.sqrt(np.mean(Bm**2))), "B_p10_over_rms": float(np.percentile(Bm, 10) / np.sqrt(np.mean(Bm**2)))}
    # particles
    pos = pool[g.wire_distance_at(pool) > 0.03][:N_PART]
    ens = monoenergetic_ensemble(DEUTERON, N_PART, 15e3, pos, np.random.default_rng(12345))
    omega_ref = E_CHARGE / DEUTERON.mass_kg * max(brms, 1e-12)
    dt = min(pr.dt_s, 0.15 / omega_ref)
    t_max = TRANSITS * 0.82 / v_D
    sample_every = max(1, int(t_max / dt / 200))
    res = run_orbits(ens, field, t_max, dt, 0.82, sample_every=sample_every)
    lost = res.loss_times
    S = {f"S{c}": float(1 - np.sum(lost <= c * 0.82 / v_D) / N_PART) for c in (4, 8, 20)}
    # core concentration: time-averaged fraction of live particles inside r<0.35 (relative to uniform shell expectation)
    core = []
    for t, x, alive in res.extra["samples"]:
        if alive.sum() > 0:
            r = np.linalg.norm(x[alive], axis=1)
            core.append(float(np.mean(r < 0.35)))
    scaled = w.scaled(1.0).with_current(pr.current_A).with_wire_radius(a)
    med = float(np.median(lost[np.isfinite(lost)])) * v_D / 0.82 if np.isfinite(lost).any() else float("inf")
    results[w.name] = {
        "meta": w.meta, "clearance": w.clearance_m(), "active_length": w.active_length_m, "circuits": len(getattr(w, "windings", [w])),
        "topology": {k: v for k, v in topo.items() if not k.startswith("per_seed")}, "field": bstats,
        "reactor": {"survival": S, "median_loss_transits": med, "core_fraction": float(np.mean(core)) if core else float("nan"),
                    "core_fraction_initial": float(np.mean(np.linalg.norm(pos, axis=1) < 0.35)),
                    "energy_drift_rel_max": res.energy_drift_rel_max, "mu_variation_median": res.mu_variation_median,
                    "counts": res.counts, "joule_power_W": scaled.joule_power_W(), "wire_radius_m": a, "dt": dt, "n": N_PART},
        "grid_n": GRID_N, "grid_seconds": t_grid, "seconds": time.time() - t0,
    }
    print(f"{w.name:52s} clr={w.clearance_m():.4f} closed={topo['frac_closed']:.2f} Rmir p90={topo['mirror_ratio_p90']:.2f} pred={topo['predicted_adiabatic_retention']:.3f} | "
          f"S4={S['S4']:.2f} S8={S['S8']:.2f} S20={S['S20']:.2f} core={np.mean(core):.2f} B_rms={brms:.2f} T Bmax/rms={bstats['B_max_over_rms']:.1f} ({time.time()-t0:.0f}s)", flush=True)
    json.dump(results, open(out_path, "w"), default=float)
print("done", part)
