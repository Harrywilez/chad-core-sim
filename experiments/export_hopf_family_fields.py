"""Field character of the exported geometries (adds a "fields" section to exports/hopf_family/MANIFEST.json).

For each geometry at 1 m scale and 1 A per strand: the toroidal field at the tube centre (ρ = R0, z = 0,
averaged over φ), the vertical field at the origin and on the midplane inside the hole, the field just
outside the torus, and the magnetic dipole moment of the whole circuit — the numbers that separate a
screw pinch (single strand), a toroidal-field coil (mirror pair, opposed currents) and a ring current
(mirror pair, equal currents), and that show what the open ends of version 1 leave behind.

    python3 experiments/export_hopf_family_fields.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ccsim import geometry as G  # noqa: E402
from ccsim.fields import biot_savart  # noqa: E402

OUT = ROOT / "exports" / "hopf_family"
M = json.load(open(OUT / "MANIFEST.json"))


def winding(name, pts, current, a, closed):
    if closed:
        return G.Winding(name, pts, current, a, closed=True, active_count=len(pts))
    P, n_active = G.close_with_return(pts, 3.8)
    return G.Winding(name, P, current, a, closed=True, active_count=n_active)


def dipole_moment(w: G.Winding) -> np.ndarray:
    p1, dl = w.segments()
    mid = p1 + 0.5 * dl
    return 0.5 * w.current_A * np.cross(mid, dl).sum(axis=0)


def probe(ws: list, R0: float, r: float) -> dict:
    """ws: list of Windings (1 m scale).  Returns field numbers in tesla per ampere."""
    phis = np.linspace(0, 2 * np.pi, 48, endpoint=False)
    centre_pts = np.column_stack((R0 * np.cos(phis), R0 * np.sin(phis), np.zeros_like(phis)))
    B = sum(biot_savart(centre_pts, w) for w in ws)
    phihat = np.column_stack((-np.sin(phis), np.cos(phis), np.zeros_like(phis)))
    Bphi = np.sum(B * phihat, axis=1)
    Bz_centre = np.sum(B[:, 2])
    mid_hole = np.array([[0.0, 0.0, 0.0], [0.5 * (R0 - r), 0.0, 0.0], [0.0, 0.0, 0.5], [0.0, 0.0, 1.0]])
    Bh = sum(biot_savart(mid_hole, w) for w in ws)
    outside = np.column_stack(((R0 + r + 0.1) * np.cos(phis), (R0 + r + 0.1) * np.sin(phis), np.zeros_like(phis)))
    Bo = sum(biot_savart(outside, w) for w in ws)
    # rotational transform proxy: the poloidal component at the tube centre relative to the toroidal one, φ-averaged
    m = sum(dipole_moment(w) for w in ws)
    return {"B_tor_tube_centre_mean": float(Bphi.mean()), "B_tor_tube_centre_ripple": float((Bphi.max() - Bphi.min()) / abs(Bphi.mean())),
            "Bz_tube_centre_mean": float(np.mean(B[:, 2])), "B_origin": Bh[0].tolist(), "B_hole_midplane_half_radius": Bh[1].tolist(),
            "B_axis_z0.5": Bh[2].tolist(), "B_axis_z1.0": Bh[3].tolist(),
            "Bmag_outside_0.1_beyond_edge_mean": float(np.linalg.norm(Bo, axis=1).mean()),
            "dipole_moment_A_m2": m.tolist(), "dipole_moment_mag": float(np.linalg.norm(m))}


def main() -> None:
    fields: dict = {"units": "tesla per ampere of strand current, 1 m unit ball; dipole moment in A·m² for 1 A; open strands closed by the evaluator's remote return (radius 3.8 m)"}
    g = M["geometries"]
    # 1. drift winding (open, with return) and the hopf torus strand
    for key, gen in (("hopf_drift_s64_180", lambda: G.hopf_drift("s64_180", 10)), ("hopf_torus_eta070_24x1", lambda: G.hopf_torus(0.70, 24, 1.0))):
        P = gen()
        fr = M["facts"]["torus_frames"]["0.70"] if "torus" in key else None
        if fr is None:
            s0 = g[key]["scale_s0"]
            R0 = 0.5 * (g[key]["torus_at_eta0"]["R0"] + g[key]["torus_at_eta1"]["R0"])
            r = 0.5 * (g[key]["torus_at_eta0"]["r"] + g[key]["torus_at_eta1"]["r"])
        else:
            R0, r = fr["R0"], fr["r"]
        w = winding(key, P, 1.0, 0.01, closed=False)
        res = probe([w], R0, r)
        res["probe_torus_R0_r"] = [R0, r]
        fields[key] = res
        print(key, {k: (round(v, 6) if isinstance(v, float) else v) for k, v in res.items() if k in ("B_tor_tube_centre_mean", "Bz_tube_centre_mean", "B_origin", "dipole_moment_mag")})
    # 2. mirror pairs: opposed and equal currents
    A, B = G.hopf_mirror_pair(0.70, 24, 1.0, "woven", 0.025, tau0=np.pi / 2, mirror_rotation=0.0)
    (A2, B2), _ = G.hopf_helical_pair(0.60, 24, eps=0.0, periods=1, delta_eta=0.015)
    (A3, B3), _ = G.hopf_helical_pair(0.60, 24, eps=0.70, periods=4, delta_eta=0.015)
    for key, (PA, PB), closed, eta, a in (("mirror_pair_v1_woven", (A, B), False, 0.70, 0.0045),
                                          ("meshed_torus_II_achiral_de015", (A2, B2), True, 0.60, 0.00448),
                                          ("meshed_torus_II_chiral_e070_n4_de015", (A3, B3), True, 0.60, 0.00448)):
        fr = M["facts"]["torus_frames"][f"{eta:.2f}"]
        for sense, lab in ((-1, "opposed"), (1, "equal")):
            ws = [winding(key + " A", PA, 1.0, a, closed), winding(key + " B", PB, float(sense), a, closed)]
            res = probe(ws, fr["R0"], fr["r"])
            if not closed:
                res["note"] = "open strands closed by the evaluator's remote return leads; the dipole moment includes the leads"
                # the leads: how close do they come to the tube's centre circle?
                ring = np.column_stack((fr["R0"] * np.cos(np.linspace(0, 2 * np.pi, 720)), fr["R0"] * np.sin(np.linspace(0, 2 * np.pi, 720)), np.zeros(720)))
                from scipy.spatial import cKDTree
                res["return_lead_min_distance_to_tube_axis"] = float(min(cKDTree(ring).query(w.points[w.active_count:])[0].min() for w in ws))
                res["tube_radius"] = fr["r"]
                ws_open = [G.Winding(key + " A", PA, 1.0, a), G.Winding(key + " B", PB, float(sense), a)]
                res_open = probe(ws_open, fr["R0"], fr["r"])
                res["open_filaments_no_leads"] = {k: res_open[k] for k in ("B_tor_tube_centre_mean", "B_tor_tube_centre_ripple", "B_origin", "Bmag_outside_0.1_beyond_edge_mean")}
            fields[f"{key}__{lab}"] = res
            print(key, lab, {k: (f"{v:.4e}" if isinstance(v, float) else [f"{x:.3e}" for x in v]) for k, v in res.items() if k in ("B_tor_tube_centre_mean", "B_tor_tube_centre_ripple", "Bz_tube_centre_mean", "B_origin", "Bmag_outside_0.1_beyond_edge_mean", "dipole_moment_mag", "return_lead_min_distance_to_tube_axis")})
    M["fields"] = fields
    json.dump(M, open(OUT / "MANIFEST.json", "w"), indent=1, default=float)


if __name__ == "__main__":
    main()
