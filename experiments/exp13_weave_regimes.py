"""exp13 — rung 1 revisited with the corrected weave: how the flux surfaces of the chiral meshed
torus depend on the over/under offset δη.

Session 8 found that the crossing-based weave let the two strands touch where the modulated
turns run nearly parallel, and that the old clearance estimate (every 3rd point) hid it.  The
weave is now built from *contact zones* (ccsim.geometry._zone_weave): over/under is frozen
wherever the strands are within h = 2 δη |∂p/∂η| of each other in the surface, so at δη ≳ 0.03
the pair is simply nested (outer strand A, inner strand B) and at δη ≲ 0.015 it is woven.  A
nested pair no longer cancels the toroidal currents locally — the gap between the shells carries
a coaxial poloidal field — so the field, the surfaces and the transform all depend on δη.  This
scan measures that dependence for the headline configurations, on the exact field, with the
same Poincaré analysis as exp10 (analyse_from_axis, 24 seeds, 40 toroidal turns).

Usage: python3 experiments/exp13_weave_regimes.py [part]
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ccsim.geometry import MultiWinding, Winding, hopf_helical_pair, hopf_torus_frame  # noqa: E402
from ccsim.surfaces import analyse_from_axis  # noqa: E402

R = ROOT / "results"
part = sys.argv[1] if len(sys.argv) > 1 else "all"


def chiral_winding(name, eta, eps, n, deta, circuits=24):
    (A, B), info = hopf_helical_pair(eta, circuits, eps=eps, periods=n, delta_eta=deta)
    ws = [Winding(f"{name} #{k}", pts, sgn, 0.005, closed=True, active_count=len(pts), meta={"index": k, "sign": sgn})
          for k, (pts, sgn) in enumerate(((A, 1), (B, -1)))]
    return MultiWinding(name, ws, {"family": "hopf-mirror-chiral", "eta": eta, "circuits": circuits, "eps": eps, "periods": n,
                                   "delta_eta": deta, "sense": -1, **info})


cases = []
for deta in (0.005, 0.01, 0.015, 0.02, 0.03, 0.05):
    cases.append((f"η 0.60, ε 0.70, n 4, δη {deta}", 0.60, 0.70, 4, deta))
for deta in (0.01, 0.03):
    cases.append((f"η 0.50, ε 0.55, n 3, δη {deta}", 0.50, 0.55, 3, deta))
    cases.append((f"η 0.60, ε 0 (pure pair), δη {deta}", 0.60, 0.0, 1, deta))
if part != "all":
    cases = [c for i, c in enumerate(cases) if i % 2 == int(part)]

out_path = R / f"exp13_{part}.json"
results = json.load(open(out_path)) if out_path.exists() else {}
for name, eta, eps, n, deta in cases:
    if name in results:
        continue
    t0 = time.time()
    w = chiral_winding(name, eta, eps, n, deta)
    R0, rt = hopf_torus_frame(eta)
    res = analyse_from_axis(w, R0, rt, n_seeds=24, max_turns=40)
    res["meta"] = w.meta
    res["clearance"] = w.clearance_m()
    res["seconds"] = time.time() - t0
    results[name] = res
    prof = " ".join(f"{r:.2f}:{i:+.3f}" for r, i in res["iota_profile"][::3])
    print(f"{name:36s} weave={w.meta['weave']:6s} zones={w.meta['zone_fraction']:.2f} clr={res['clearance']:.4f} surf={res['frac_surface']:.2f} "
          f"open={res['frac_open']:.2f} r_out={res['outermost_surface_r']:.3f} ι_axis={res['iota_axis']:+.3f} ι_edge={res['iota_edge']:+.3f} "
          f"well={res['well_depth']:+.3f} | {prof} ({time.time()-t0:.0f}s)", flush=True)
    json.dump(results, open(out_path, "w"), default=float)
print("done", part)
