"""Field-line topology of every configuration and the adiabatic confinement prediction.

Seeds are the same 256 initial positions the particle screen uses (unit ball).
Outputs per winding: fraction of lines ending on the wall / a conductor / closed,
connection length, mirror ratio distribution, predicted isotropic trapped
fraction sqrt(1 - B_seed/B_mirror).  The prediction is tested against the
reactor-scale D+ retention (adiabaticity ~0.01) in exp04.
"""
import json, sys, time
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ccsim.fieldlines import analyse, summarise
from ccsim.runner import build_grids, initial_positions, unit_catalogue

cat = unit_catalogue()
grids = build_grids(cat, ROOT / "cache", verbose=False)
rng = np.random.default_rng(271828)
union = np.vstack([w.active_points[::4] for w in cat])
seeds = initial_positions(256, rng, union)
out = {}
for w in cat:
    t0 = time.time()
    lines = analyse(grids[w.name], seeds[:128], wall_radius=0.82, step=0.006, max_length=40.0, wire_clear=0.008)
    s = summarise(lines)
    s["per_seed_predicted"] = [l.predicted_trapped_fraction for l in lines]
    s["per_seed_mirror_ratio"] = [l.mirror_ratio if np.isfinite(l.mirror_ratio) else 1e9 for l in lines]
    s["per_seed_ends"] = [[l.end_plus, l.end_minus] for l in lines]
    out[w.name] = s
    print(f"{w.name:55s} wall-both={s['frac_wall_both']:.2f} wire={s['frac_wire_any']:.2f} closed={s['frac_closed']:.2f} "
          f"L_conn={s['connection_length_median']:.2f} R_mirror med={s['mirror_ratio_median']:.2f} p90={s['mirror_ratio_p90']:.2f} "
          f"pred adiabatic retention={s['predicted_adiabatic_retention']:.3f} ({time.time()-t0:.0f}s)", flush=True)
json.dump(out, open(ROOT / "results" / "exp03_topology.json", "w"))
