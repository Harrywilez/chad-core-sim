"""Cheap pre-screen: geometry (clearance / r_max / points) and, with --field, the
Biot-Savart grid only (no particles, no score) so I can predict the buildability
factor before spending an evaluation.  Reuses the same grid cache as ccsim.evaluate,
so a pre-screened design costs nothing extra when it is later evaluated.

Usage:
  python3 agent/prescreen.py designs/x.json [designs/y.json ...] [--field]
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ccsim.design import build, validate  # noqa: E402
from ccsim.fields import FieldGrid  # noqa: E402
from ccsim.forces import CONDUCTOR_LIMITS_A_M2  # noqa: E402
from ccsim.presets import PRESETS  # noqa: E402
from ccsim.runner import make_field, roi_mask  # noqa: E402

CACHE = ROOT / "cache" / "designs"
J_LIM = CONDUCTOR_LIMITS_A_M2["REBCO tape at 20 T, 20 K (engineering)"]


def grid_for(design, grid_n=25):
    v = validate(design)
    w = build(design)
    CACHE.mkdir(parents=True, exist_ok=True)
    gpath = CACHE / f"grid_{v['hash']}_n{grid_n}.npz"
    if gpath.exists():
        d = np.load(gpath)
        g = FieldGrid.__new__(FieldGrid)
        g.winding, g.half, g.n = w, 0.87, grid_n
        g.axis = np.linspace(-0.87, 0.87, grid_n)
        g.dx = float(g.axis[1] - g.axis[0])
        xx, yy, zz = np.meshgrid(g.axis, g.axis, g.axis, indexing="ij")
        g.points = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
        g.B = d["B"]
        g.A = None
        g.Bmag = np.linalg.norm(g.B, axis=3)
        g.wire_distance = d["wire_distance"]
    else:
        g = FieldGrid(w, 0.87, grid_n, with_A=False)
        np.savez_compressed(gpath, B=g.B, wire_distance=g.wire_distance)
    return v, w, g


def report(design, with_field=False, B_target=2.0):
    v = validate(design)
    out = {"name": design.get("name"), "valid": v["valid"], "clearance": v["clearance"], "r_max": v["r_max"],
           "points": v["conductor_points"], "circuits": v["circuits"], "problems": v["problems"],
           "wr_best": v["clearance"] / 3.0}
    if with_field and v["valid"]:
        v, w, g = grid_for(design)
        pr = PRESETS["reactor"]()
        mask = roi_mask(g)
        Bm = g.Bmag.ravel()[mask]
        field0, a = make_field(g, pr, w, pr.current_A, pr.wire_radius_m)
        brms0 = float(abs(field0.gain) * np.sqrt(np.mean(Bm ** 2)))
        current = pr.current_A * B_target / max(brms0, 1e-12)
        J = current / (math.pi * a * a)
        wr = float(design.get("wire_radius", 0.005))
        b_clear = min(1.0, v["clearance"] / (3 * wr))
        b_j = min(1.0, J_LIM / J)
        out.update({"a_used": a, "current_kA": current / 1e3, "J_over_REBCO": J / J_LIM,
                    "b_clear": b_clear, "b_J": b_j, "build": b_clear * b_j,
                    "build_max": min(1.0, J_LIM / J),  # with wr set to clearance/3
                    "Bmax_over_rms": float(Bm.max() / np.sqrt(np.mean(Bm ** 2)))})
    return out


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    fld = "--field" in sys.argv
    for p in args:
        d = json.load(open(p))
        r = report(d, with_field=fld)
        print(f"{Path(p).name:34s} clr={r['clearance']:.4f} rmax={r['r_max']:.3f} pts={r['points']:6d} "
              f"wr*={r['wr_best']:.4f} valid={r['valid']}"
              + (f" I={r.get('current_kA', float('nan')):7.1f}kA J/Jr={r.get('J_over_REBCO', float('nan')):6.3f} "
                 f"build={r.get('build', float('nan')):.3f} build_max={r.get('build_max', float('nan')):.3f}" if fld else "")
              + ("" if r["valid"] else f"  {r['problems']}"))
