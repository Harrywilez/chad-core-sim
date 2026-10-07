"""exp11 — rung 2: guiding-centre confinement inside the flux-surface region, with and
without collisions, untwisted vs chiral (stellarator) meshed torus.

Particles: 15 keV D⁺ at the reactor preset (1 m ball, 364 kA per strand), seeded inside
0.7 × (last closed surface radius) around the magnetic axis found by rung 1, isotropic
pitch.  Loss = reaching the sphere r = 1.05 or coming within 0.02 of a strand (the
winding is the vessel).  Guiding-centre RK4 on the exact Biot–Savart field.
  ν = 0                   collisionless: drift losses (untwisted) / ripple-trapped losses (twisted)
  ν = 0.1, 0.5 / transit  pitch-angle scattering scaled up ~10³–10⁴ × the physical rate,
                          to see the collisional (1/ν) channel inside a tractable run

Usage: python3 experiments/exp11_transport.py <run index 0..5 | all>
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ccsim.constants import DEUTERON, E_CHARGE
from ccsim.fastfield import SegmentSource
from ccsim.geometry import MultiWinding, Winding, hopf_helical_pair, hopf_mirror_pair, hopf_torus_frame
from ccsim.guiding_center import run_gc
from ccsim.presets import PRESETS
from ccsim.surfaces import find_axis

R = ROOT / "results"
which = sys.argv[1] if len(sys.argv) > 1 else "all"
pr = PRESETS["reactor"]()


def pair_to_winding(name, A, B, meta):
    ws = [Winding(f"{name} #{k}", pts, sgn, 0.005, closed=True, active_count=len(pts)) for k, (pts, sgn) in enumerate(((A, 1), (B, -1)))]
    return MultiWinding(name, ws, meta)


def build(kind, eta, eps, n):
    if kind == "untwisted":
        A, B = hopf_mirror_pair(eta, 24, revolutions=2.0, delta_eta=0.03)
        info = {}
    else:
        (A, B), info = hopf_helical_pair(eta, 24, eps=eps, periods=n, delta_eta=0.03, helical_mode="toroidal")
    return pair_to_winding(f"{kind} η={eta} ε={eps} n={n}", A, B, {"kind": kind, "eta": eta, "eps": eps, "periods": n, **info})


runs = [
    ("untwisted η=0.60", "untwisted", 0.60, 0.0, 4, 0.0, 0.11),
    ("chiral ε=0.70 n=4 η=0.60", "chiral", 0.60, 0.70, 4, 0.0, 0.11),
    ("chiral ε=0.55 n=3 η=0.50", "chiral", 0.50, 0.55, 3, 0.0, 0.085),
    ("chiral ε=0.70 n=4 η=0.60, ν=0.1/transit", "chiral", 0.60, 0.70, 4, 0.1, 0.11),
    ("chiral ε=0.70 n=4 η=0.60, ν=0.5/transit", "chiral", 0.60, 0.70, 4, 0.5, 0.11),
    ("untwisted η=0.60, ν=0.1/transit", "untwisted", 0.60, 0.0, 4, 0.1, 0.11),
]
sel = runs if which == "all" else [runs[int(which)]]
out_path = R / f"exp11_{which}.json"
results = json.load(open(out_path)) if out_path.exists() else {}
N_PART, TRANSITS = 64, 60
v = math.sqrt(2 * 15e3 * E_CHARGE / DEUTERON.mass_kg)
transit = 0.82 / v

for name, kind, eta, eps, n, nu_per_transit, r_seed in sel:
    if name in results:
        continue
    t0 = time.time()
    w = build(kind, eta, eps, n)
    clr = w.clearance_m()
    a = min(pr.wire_radius_m, 0.45 * clr)
    ws = w.with_current(pr.current_A).with_wire_radius(a)
    src = SegmentSource(ws)
    R0, rt = hopf_torus_frame(eta)
    axis = np.array([R0, 0.0]) if kind == "untwisted" else find_axis(ws, R0, rt)
    rng = np.random.default_rng(7)
    th = rng.uniform(0, 2 * np.pi, N_PART); ze = rng.uniform(0, 2 * np.pi, N_PART); rr = r_seed * np.sqrt(rng.uniform(0, 1, N_PART))
    # seed disc around the axis at φ = ζ, rotated toroidally (the axis position at other ζ differs slightly for the chiral case; fine for a seed)
    X0 = np.column_stack(((axis[0] + rr * np.cos(th)) * np.cos(ze), (axis[0] + rr * np.cos(th)) * np.sin(ze), axis[1] + rr * np.sin(th)))
    d = rng.normal(size=(N_PART, 3)); d /= np.linalg.norm(d, axis=1)[:, None]
    tree = cKDTree(ws.active_points[::2])

    def loss_fn(X):
        return (np.linalg.norm(X, axis=1) > 1.05) | (tree.query(X)[0] < 0.02)

    res = run_gc(src, X0, v * d, DEUTERON, TRANSITS * transit, 0.004 / v, loss_fn, nu=nu_per_transit / transit, rng=rng, sample_every=25)
    lt = res.loss_times
    S = {f"S{c}": float(np.mean(lt > c * transit)) for c in (5, 10, 20, 40, 60)}
    results[name] = {"meta": w.meta, "axis": axis.tolist(), "r_seed": r_seed, "nu_per_transit": nu_per_transit, "n": N_PART, "transits": TRANSITS,
                     "survival": S, "curve_t_transits": (res.t / transit).tolist(), "curve_alive": res.alive_frac.tolist(),
                     "loss_times_transits": np.where(np.isfinite(lt), lt / transit, np.inf).tolist(), "energy_drift": res.energy_drift_rel_max,
                     "clearance": clr, "wire_radius_m": a, "seconds": time.time() - t0}
    print(f"{name:45s} " + " ".join(f"{k}={v_:.2f}" for k, v_ in S.items()) + f"  axis=({axis[0]:.3f},{axis[1]:.3f}) ({time.time()-t0:.0f}s)", flush=True)
    json.dump(results, open(out_path, "w"), default=float)
print("done", which)
