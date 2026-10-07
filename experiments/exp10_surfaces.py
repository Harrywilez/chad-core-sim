"""exp10 — rung 1: does the meshed torus have flux surfaces, and can a twist of the weave
give it rotational transform?

Every configuration is the woven Hopf mirror pair (η 0.70, 24 circuits, opposed currents,
strand phase tau0 = 0) on a winding surface deformed by torus_deformation:
ellipticity ε (rotating ellipse, l = 2), helical-axis excursion h, n field periods.
Field lines are traced on the exact Biot–Savart field (numba) from 16 seeds on the
outboard ray at φ = 0, 30 toroidal turns each; the loss boundary is the winding itself
(2 % of the ball from any strand) — the vessel of a toroidal device is inside the coil,
not the sphere used for the single-particle screens.

Usage: python3 experiments/exp10_surfaces.py <part 0|1|all>
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ccsim.geometry import MultiWinding, Winding, close_with_return, hopf_torus_frame, twisted_mirror_pair
from ccsim.surfaces import analyse_surfaces, seeds_on_ray

R = ROOT / "results"
part = sys.argv[1] if len(sys.argv) > 1 else "all"


def pair_winding(name, eps, n, h, tri=0.0, deta=None, circuits=24, eta=0.70):
    if deta is None:
        deta = 0.025 + 0.08 * (eps + h + tri)
    A, B = twisted_mirror_pair(eta, circuits, ellipticity=eps, periods=n, helical_axis=h, triangularity=tri, delta_eta=deta)
    ws = []
    for k, (pts, sgn) in enumerate(((A, 1), (B, -1))):
        closed, na = close_with_return(pts, 3.8)
        ws.append(Winding(f"{name} #{k}", closed, sgn, 0.005, closed=True, active_count=na, meta={"index": k, "sign": sgn}))
    return MultiWinding(name, ws, {"family": "hopf-mirror-twisted", "eta": eta, "circuits": circuits, "ellipticity": eps, "periods": n,
                                   "helical_axis": h, "triangularity": tri, "delta_eta": deta, "sense": -1})


def chiral_winding(name, eps, n, deta=0.03, circuits=24, eta=0.70, deform_eps=0.0, helical_mode="toroidal", helical_axis=0.0):
    from ccsim.geometry import hopf_helical_pair, torus_deformation, hopf_torus_frame, _hopf_points
    deform = None
    if deform_eps or helical_axis:
        R0_, rt_ = hopf_torus_frame(eta)
        u_ = np.linspace(0, 1, 4001); raw = _hopf_points(np.full_like(u_, eta), np.zeros_like(u_), 2 * np.pi * u_); s_ = np.linalg.norm(raw, axis=1).max()
        deform = torus_deformation(R0_ * s_, rt_ * s_, deform_eps, n, helical_axis)
    (A, B), info = hopf_helical_pair(eta, circuits, eps=eps, periods=n, delta_eta=deta, deform=deform, helical_mode=helical_mode)
    ws = [Winding(f"{name} #{k}", pts, sgn, 0.005, closed=True, active_count=len(pts), meta={"index": k, "sign": sgn}) for k, (pts, sgn) in enumerate(((A, 1), (B, -1)))]
    return MultiWinding(name, ws, {"family": "hopf-mirror-chiral", "eta": eta, "circuits": circuits, "eps": eps, "periods": n, "delta_eta": deta,
                                   "deform_eps": deform_eps, "helical_mode": helical_mode, "sense": -1, **info})


cases = [
    ("baseline (round torus)", 0.0, 5, 0.0, 0.0),
    ("ellipse 0.10, n=5", 0.10, 5, 0.0, 0.0),
    ("ellipse 0.20, n=5", 0.20, 5, 0.0, 0.0),
    ("ellipse 0.30, n=5", 0.30, 5, 0.0, 0.0),
    ("ellipse 0.20, n=3", 0.20, 3, 0.0, 0.0),
    ("ellipse 0.20, n=7", 0.20, 7, 0.0, 0.0),
    ("helical axis 0.15, n=5", 0.0, 5, 0.15, 0.0),
    ("helical axis 0.30, n=5", 0.0, 5, 0.30, 0.0),
    ("ellipse 0.20 + helical axis 0.15, n=5", 0.20, 5, 0.15, 0.0),
    ("ellipse 0.20 + helical axis 0.30, n=5", 0.20, 5, 0.30, 0.0),
    ("ellipse 0.20 + triangle 0.15, n=5", 0.20, 5, 0.0, 0.15),
    ("ellipse 0.30 + helical axis 0.15, n=4", 0.30, 4, 0.15, 0.0),
]
chiral_cases = [("chiral ε=0.20, n=5", 0.20, 5, 0.0), ("chiral ε=0.40, n=5", 0.40, 5, 0.0), ("chiral ε=0.20, n=3", 0.20, 3, 0.0),
                ("chiral ε=0.40, n=3", 0.40, 3, 0.0), ("chiral ε=0.60, n=5", 0.60, 5, 0.0), ("chiral ε=0.40, n=5 + ellipse 0.15", 0.40, 5, 0.15),
                ("chiral ε=0.40, n=7", 0.40, 7, 0.0), ("chiral ε=0.80, n=5", 0.80, 5, 0.0)]
helical_cases = [("helical-current ε=0.10, n=5", 0.10, 5), ("helical-current ε=0.20, n=5", 0.20, 5), ("helical-current ε=0.30, n=5", 0.30, 5),
                 ("helical-current ε=0.20, n=3", 0.20, 3), ("helical-current ε=0.20, n=7", 0.20, 7), ("helical-current ε=0.40, n=5", 0.40, 5),
                 ("helical-current ε=0.30, n=4", 0.30, 4), ("helical-current ε=0.30, n=3", 0.30, 3),
                 ("helical-current ε=0.55, n=5", 0.55, 5), ("helical-current ε=0.70, n=5", 0.70, 5), ("helical-current ε=0.55, n=4", 0.55, 4),
                 ("helical-current ε=0.70, n=4", 0.70, 4), ("helical-current ε=0.55, n=6", 0.55, 6), ("helical-current ε=0.85, n=5", 0.85, 5)]
best_cases = [("best: ε=0.70, n=4", 0.70, 4, 0.0, 0.70), ("best: ε=0.70, n=3", 0.70, 3, 0.0, 0.70), ("best: ε=0.70, n=4, η=0.60", 0.70, 4, 0.0, 0.60),
              ("best: ε=0.70, n=4, η=0.80", 0.70, 4, 0.0, 0.80), ("best: ε=0.55, n=4, η=0.60", 0.55, 4, 0.0, 0.60),
              ("best: ε=0.70, n=5, η=0.60", 0.70, 5, 0.0, 0.60), ("best: ε=0.85, n=4", 0.85, 4, 0.0, 0.70), ("best: ε=0.70, n=4, 36 circuits", 0.70, 4, 0.0, 0.70),
              ("best: ε=0.70, n=4, η=0.50", 0.70, 4, 0.0, 0.50), ("best: ε=0.70, n=3, η=0.60", 0.70, 3, 0.0, 0.60), ("best: ε=0.55, n=3, η=0.50", 0.55, 3, 0.0, 0.50),
              ("best: ε=0.40, n=4, η=0.60", 0.40, 4, 0.0, 0.60), ("best: baseline η=0.60 (ε=0)", 0.0, 4, 0.0, 0.60), ("best: ε=0.85, n=4, η=0.60", 0.85, 4, 0.0, 0.60)]
if part.startswith("best"):
    k = int(part[-1]) if part[-1].isdigit() else None
    cases = [(n, e, per, h, 0.0, 0.0, "toroidal", eta_) for n, e, per, h, eta_ in best_cases]
    if k is not None:
        cases = [c for i, c in enumerate(cases) if i % 2 == k]
elif part.startswith("helical"):
    k = int(part[-1]) if part[-1].isdigit() else None
    cases = [(n, e, per, 0.0, 0.0, 0.0, "toroidal") for n, e, per in helical_cases]
    if k is not None:
        cases = [c for i, c in enumerate(cases) if i % 2 == k]
elif part.startswith("chiral"):
    k = int(part[-1]) if part[-1].isdigit() else None
    cases = [(n, e, per, 0.0, 0.0, de, "poloidal") for n, e, per, de in chiral_cases]
    if k is not None:
        cases = [c for i, c in enumerate(cases) if i % 2 == k]
elif part != "all":
    cases = [c for i, c in enumerate(cases) if i % 2 == int(part)]

out_path = R / f"exp10_{part}.json"
results = json.load(open(out_path)) if out_path.exists() else {}
R0, rt = hopf_torus_frame(0.70)
for case in cases:
    name, eps, n, h, tri = case[:5]
    if name in results:
        continue
    t0 = time.time()
    chiral = len(case) > 5
    eta_c = case[7] if len(case) > 7 else 0.70
    circuits_c = 36 if "36 circuits" in name else 24
    if chiral:
        w = chiral_winding(name, eps, n, deform_eps=case[5], helical_mode=case[6], helical_axis=h, eta=eta_c, circuits=circuits_c)
        R0c, rtc = hopf_torus_frame(eta_c)
    else:
        w = pair_winding(name, eps, n, h, tri)
        R0c, rtc = R0, rt
    axis0 = np.array([R0c + h * rtc, 0.0])
    if part.startswith("best"):
        from ccsim.surfaces import analyse_from_axis
        res = analyse_from_axis(w, R0c, rtc, n_seeds=24, max_turns=40)
    else:
        seeds = seeds_on_ray(axis0[0], rtc * (1 - (0.0 if chiral else eps)), 16, s_min=0.03, s_max=0.85)
        res = analyse_surfaces(w, seeds, axis0, step=0.006, max_turns=30, wall=1.05, wire_clear=0.02)
    res["meta"] = w.meta
    res["clearance"] = w.clearance_m()
    res["seconds"] = time.time() - t0
    results[name] = res
    prof = " ".join(f"{r:.2f}:{i:+.3f}" for r, i in res["iota_profile"][::3])
    print(f"{name:42s} clr={res['clearance']:.4f} surf={res['frac_surface']:.2f} isl={res['frac_island_chaotic']:.2f} open={res['frac_open']:.2f} "
          f"r_out={res['outermost_surface_r']:.3f} ι_axis={res['iota_axis']:+.3f} ι_edge={res['iota_edge']:+.3f} well={res['well_depth']:+.3f} | {prof} ({time.time()-t0:.0f}s)", flush=True)
    json.dump(results, open(out_path, "w"), default=float)
print("done", part)
