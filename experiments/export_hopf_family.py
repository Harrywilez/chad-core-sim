"""Pull out the Hopf-drift winding and the two versions of the Hopf torus + mirror image as files:
centreline polylines (OBJ, CSV), swept tube meshes (binary STL, OBJ+MTL with one coloured object per
strand), crossing tables, and a MANIFEST.json with every parameter and the geometric numbers quoted in
exports/hopf_family/README.md.

    python3 experiments/export_hopf_family.py            # everything (~1 min)
    python3 experiments/export_hopf_family.py --no-tubes # centrelines + numbers only

Coordinates are the generators' unit-ball units (multiply by the vessel scale, e.g. 1 m).
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ccsim import geometry as G  # noqa: E402
from ccsim.export3d import (mesh_is_closed, tube_mesh, write_csv, write_obj_mesh, write_obj_polylines,  # noqa: E402
                            write_stl)

OUT = ROOT / "exports" / "hopf_family"
COL_A = (0.16, 0.45, 0.85)      # strand A blue
COL_B = (0.86, 0.30, 0.20)      # strand B (mirror) red
COL_DRIFT = (0.20, 0.55, 0.35)


def frame(eta: float) -> dict:
    """Round torus that the Hopf torus η projects to: raw (R0 = sec η, r = tan η) and normalised so the
    outer equator sits at radius 1 (R0 = 1/(1 + sin η), r = sin η/(1 + sin η))."""
    s = math.sin(eta)
    return {"R0_raw": 1 / math.cos(eta), "r_raw": math.tan(eta), "R0": 1 / (1 + s), "r": s / (1 + s),
            "hole_radius": (1 - s) / (1 + s), "aspect_ratio_R0_over_r": 1 / s}


def fit_circle(pts: np.ndarray) -> dict:
    """Least-squares circle through (near-)planar points: plane normal, centre, radius, rms residual."""
    c0 = pts.mean(axis=0)
    _, _, Vt = np.linalg.svd(pts - c0)
    n = Vt[2]
    e1, e2 = Vt[0], Vt[1]
    x, y = (pts - c0) @ e1, (pts - c0) @ e2
    Amat = np.column_stack((2 * x, 2 * y, np.ones_like(x)))
    sol, *_ = np.linalg.lstsq(Amat, x * x + y * y, rcond=None)
    cx, cy, k = sol
    R = math.sqrt(k + cx * cx + cy * cy)
    centre = c0 + cx * e1 + cy * e2
    resid = np.hypot(x - cx, y - cy) - R
    return {"normal": n.tolist(), "centre": centre.tolist(), "radius": R, "rms_residual": float(np.sqrt(np.mean(resid**2))),
            "out_of_plane_max": float(np.abs((pts - c0) @ n).max())}


def strand_numbers(P: np.ndarray, closed: bool) -> dict:
    L = G.polyline_length(P)
    clr, _ = G.min_nonlocal_distance(P, skip=12)
    return {"points": int(len(P)), "length": float(L), "closed": bool(closed), "start": P[0].tolist(), "end": P[-1].tolist(),
            "end_gap": float(np.linalg.norm(P[-1] - P[0])), "max_radius": float(np.linalg.norm(P, axis=1).max()),
            "self_clearance": float(clr)}


def write_mtl(path, materials: dict) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        for name, (r, g, b) in materials.items():
            fh.write(f"newmtl {name}\nKd {r:.3f} {g:.3f} {b:.3f}\nKa {0.2*r:.3f} {0.2*g:.3f} {0.2*b:.3f}\nKs 0.25 0.25 0.25\nNs 40\nd 1.0\n\n")


def write_obj_mesh_mtl(path: Path, parts, colors, comment=""):
    """OBJ mesh with a companion .mtl so each strand gets its own colour in viewers that honour materials."""
    mtl = path.with_suffix(".mtl")
    write_mtl(mtl, {name: col for (name, _, _), col in zip(parts, colors)})
    with open(path, "w", encoding="utf-8") as fh:
        for line in comment.splitlines():
            fh.write(f"# {line}\n")
        fh.write(f"mtllib {mtl.name}\n")
        off = 1
        for name, V, F in parts:
            fh.write(f"o {name}\nusemtl {name}\n")
            fh.write("".join(f"v {x:.6f} {y:.6f} {z:.6f}\n" for x, y, z in V))
            fh.write("".join(f"f {a + off} {b + off} {c + off}\n" for a, b, c in F))
            off += len(V)


def export_strands(folder: Path, key: str, strands: list, radius: float, tubes: bool, n_seg: int, comment: str) -> dict:
    """strands: [(name, points, closed, color)].  Writes centreline OBJ + CSVs (+ tube STL/OBJ) and returns file list."""
    folder.mkdir(parents=True, exist_ok=True)
    files = {}
    write_obj_polylines(folder / f"{key}.centreline.obj", [(n, P) for n, P, _, _ in strands], comment)
    files["centreline_obj"] = f"{key}.centreline.obj"
    for n, P, closed, _ in strands:
        write_csv(folder / f"{key}.{n}.csv", P, comment + f"\nstrand {n}, {'closed' if closed else 'open'} polyline, {len(P)} points")
        files[f"csv_{n}"] = f"{key}.{n}.csv"
    if tubes:
        parts = []
        for n, P, closed, _ in strands:
            V, F = tube_mesh(P, radius, n_seg, closed=closed)
            assert mesh_is_closed(V, F), (key, n)
            parts.append((n, V, F))
        write_stl(folder / f"{key}.tube.stl", parts, header=f"{key} tube r={radius}")
        write_obj_mesh_mtl(folder / f"{key}.tube.obj", parts, [c for _, _, _, c in strands], comment + f"\ntube radius {radius}")
        files["tube_stl"] = f"{key}.tube.stl"
        files["tube_obj"] = f"{key}.tube.obj"
        files["tube_triangles"] = int(sum(len(F) for _, _, F in parts))
    return files


def crossing_points_from_params(P_fine_fn, uc):
    return P_fine_fn(np.asarray(uc))


# ---------------------------------------------------------------------------


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-tubes", action="store_true")
    a = ap.parse_args()
    tubes = not a.no_tubes
    t_start = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    manifest: dict = {"units": "generator unit-ball units (× vessel scale, 1 m in the evaluator)", "geometries": {}}

    # ======================================================================
    # 0. Hopf coordinates / fibre facts (numbers for the note)
    # ======================================================================
    facts: dict = {"torus_frames": {f"{eta:.2f}": frame(eta) for eta in (0.48, 0.60, 0.64, 0.70, 0.82)}}
    u = np.linspace(0.0, 1.0, 4001)
    eta = 0.70
    fibre = G._hopf_points(np.full_like(u, eta), np.zeros_like(u), 2 * np.pi * u)          # one Hopf fibre, φ = 0, raw frame
    anti = fibre * np.array([1.0, 1.0, -1.0])                                                # its z-mirror: an anti-fibre
    fc = fit_circle(fibre)
    fr = frame(eta)
    d = cKDTree(fibre).query(anti)[0]
    # the two families meet twice: count sign changes of the distance minima (local minima below 1e-3)
    idx = np.where((d[1:-1] < d[:-2]) & (d[1:-1] < d[2:]) & (d[1:-1] < 2e-3))[0]
    facts["fibre_eta_0.70"] = {"is_planar_circle": fc["out_of_plane_max"] < 1e-9 and fc["rms_residual"] < 1e-9,
                               "circle_radius_raw": fc["radius"], "R0_raw": fr["R0_raw"],
                               "circle_centre_distance_from_axis_raw": float(np.hypot(*fc["centre"][:2])), "r_raw": fr["r_raw"],
                               "plane_tilt_deg": math.degrees(math.acos(abs(fc["normal"][2]))), "eta_deg": math.degrees(eta),
                               "meets_its_mirror_image_n_times": int(len(idx))}
    manifest["facts"] = facts

    # ======================================================================
    # 1. Hopf drift (Codex factorial variants) + the continuation + the torus it heads for
    # ======================================================================
    f1 = OUT / "1_hopf_drift"
    for key, spec in G.HOPF_DRIFT_VARIANTS.items():
        P = G.hopf_drift(key, 10)
        Pf = G.hopf_drift(key, 10, n_per_circuit=300)
        eta0 = 0.82
        s0 = (1 - math.sin(eta0)) / math.cos(eta0)
        de = spec["eta1"] - eta0
        dphi = math.radians(spec["sweep"])
        v = np.linspace(0, 1, 4097)
        speed = np.sqrt(de * de + 0.25 * np.sin(2 * (eta0 + de * v)) ** 2 * dphi * dphi)
        base_len = float(np.sum(0.5 * (speed[1:] + speed[:-1]) * np.diff(v)))
        eta_r, a1_r, a2_r = G.hopf_coordinates(P, s0)
        num = strand_numbers(P, False)
        num.update({"variant": key, "label": spec["label"], "eta0": eta0, "eta1": spec["eta1"], "sweep_deg": spec["sweep"], "circuits": 10,
                    "tau0": "pi/2", "scale_s0": s0, "base_sphere_path_length": base_len,
                    "torus_at_eta0": {k: s0 * frame(eta0)[k + "_raw"] for k in ("R0", "r")},
                    "torus_at_eta1": {k: s0 * frame(spec["eta1"])[k + "_raw"] for k in ("R0", "r")},
                    "recovered_eta_range": [float(eta_r.min()), float(eta_r.max())],
                    "z_range": [float(P[:, 2].min()), float(P[:, 2].max())],
                    "conductor_radius_used_for_tube": 0.012})
        do_tube = tubes and key in ("s64_180", "d48_180")
        files = export_strands(f1, f"hopf_drift_{key}_10circuits", [("drift", Pf if do_tube else P, False, COL_DRIFT)] if do_tube else [("drift", P, False, COL_DRIFT)],
                               0.012, do_tube, 8, f"Hopf drift winding {key} ({spec['label']}), 10 circuits, unit ball; ccsim.geometry.hopf_drift('{key}', 10)")
        if do_tube:   # the centreline files should carry the canonical 1501-point polyline, the tube the finer one
            write_obj_polylines(f1 / f"hopf_drift_{key}_10circuits.centreline.obj", [("drift", P)], f"Hopf drift winding {key}, 10 circuits (canonical 1501-point polyline)")
            write_csv(f1 / f"hopf_drift_{key}_10circuits.drift.csv", P, f"Hopf drift winding {key}, 10 circuits (canonical 1501-point polyline)")
        num["files"] = files
        manifest["geometries"][f"hopf_drift_{key}"] = num
        print(f"drift {key}: L={num['length']:.2f} clearance={num['self_clearance']:.4f} eta {eta_r.min():.3f}..{eta_r.max():.3f} z {P[:,2].min():+.3f}..{P[:,2].max():+.3f}")
    # continuation and limit
    Pc = G.hopf_continued(0.82, 0.64, 720.0, 48)
    num = strand_numbers(Pc, False)
    num.update({"eta0": 0.82, "eta1": 0.64, "sweep_deg": 720.0, "circuits": 48, "note": "same drift law, φ sweeps 720°, normalised to the unit ball"})
    num["files"] = export_strands(f1, "hopf_continued_720deg_48circuits", [("drift", Pc, False, COL_DRIFT)], 0.008, False, 8,
                                  "Hopf drift continued: η 0.82→0.64 while φ sweeps 720°, 48 circuits; ccsim.geometry.hopf_continued()")
    manifest["geometries"]["hopf_continued"] = num
    Pt = G.hopf_torus(0.70, 24, 1.0)
    num = strand_numbers(Pt, False)
    num.update({"eta": 0.70, "circuits": 24, "revolutions": 1.0, "torus": frame(0.70), "note": "the limit of the drift: η fixed, 24 fibre turns, φ sweeps 360° (open (24.5, 23.5) strand); this is strand A of mirror pair v1 before the weave displacement"})
    num["files"] = export_strands(f1, "hopf_torus_eta070_24x1", [("torus", Pt, False, COL_A)], 0.006, False, 8,
                                  "Hopf torus strand η 0.70, 24 circuits, 1 revolution; ccsim.geometry.hopf_torus(0.70, 24, 1.0)")
    manifest["geometries"]["hopf_torus_eta070_24x1"] = num

    # ======================================================================
    # 2. Mirror pair v1 — session-5 scheme: open strands, τ0 = π/2, mirror rotation 0, (−1)^k lattice weave
    # ======================================================================
    f2 = OUT / "2_mirror_pair_v1"
    f2.mkdir(parents=True, exist_ok=True)
    eta, N, r_rev, de = 0.70, 24, 1.0, 0.025
    for mode in ("woven", "nested"):
        info: dict = {}
        A, B = G.hopf_mirror_pair(eta, N, r_rev, mode, de, tau0=np.pi / 2, mirror_rotation=0.0, info=info)
        Af, Bf = G.hopf_mirror_pair(eta, N, r_rev, mode, de, n_per_circuit=200, tau0=np.pi / 2, mirror_rotation=0.0)
        sc = info["strand_scale"]
        eA, a1A, a2A = G.hopf_coordinates(A, sc)
        eB, a1B, a2B = G.hopf_coordinates(B, sc)
        MA = A * np.array([1, 1, -1.0])
        dAB = cKDTree(A).query(B)[0].min()
        num = {"eta": eta, "circuits": N, "revolutions": r_rev, "mode": mode, "delta_eta": de, "tau0": "pi/2", "mirror_rotation": 0.0,
               "p": info["p"], "q": info["q"], "strand_scale": sc, "torus": frame(eta),
               "A": strand_numbers(A, False), "B": strand_numbers(B, False),
               "min_distance_A_B": float(dAB), "B_equals_zmirror_of_A_maxdev": float(np.abs(B - MA).max()),
               "B_equals_zmirror_of_A_nearest": float(cKDTree(MA).query(B)[0].max()),
               "eta_range_A": [float(eA.min()), float(eA.max())], "eta_range_B": [float(eB.min()), float(eB.max())],
               "conductor_radius_used_for_tube": 0.0045}
        if mode == "woven":
            uc, vc = np.array(info["crossing_u"]), np.array(info["crossing_v"])
            jj, kk, sg = np.array(info["crossing_j"]), np.array(info["crossing_k"]), np.array(info["crossing_sign"])
            # alternation along A (sorted by u) and along B (sorted by v)
            oA, oB = np.argsort(uc), np.argsort(vc)
            altA = float(np.mean(sg[oA][1:] * sg[oA][:-1] < 0))
            altB = float(np.mean(sg[oB][1:] * sg[oB][:-1] < 0))
            # positions of the crossings on the undisplaced torus (raw frame → normalised)
            p, q = info["p"], info["q"]
            X = G._hopf_points(np.full_like(uc, eta), 2 * np.pi * r_rev * uc, np.pi / 2 + 2 * np.pi * N * uc) * sc
            # check: the mirror strand's undisplaced point at v is the same place
            Y = G._hopf_points(np.full_like(vc, eta), 2 * np.pi * r_rev * vc, np.pi / 2 + 2 * np.pi * N * vc) * sc * np.array([1, 1, -1.0])
            num.update({"crossings": int(len(uc)), "two_p_q": 2 * p * q, "alternation_along_A": altA, "alternation_along_B": altB,
                        "exact_mirror_scheme_used": info["exact_mirror"], "crossing_position_check_max": float(np.abs(X - Y).max()),
                        "crossing_spacing_u": float(np.median(np.diff(np.sort(uc)))), "expected_spacing_1_over_2pq": 1 / (2 * p * q)})
            with open(f2 / f"mirror_pair_v1_{mode}.crossings.csv", "w", encoding="utf-8") as fh:
                fh.write("# crossings of strand A(u) with the mirror strand B(v): u = (j/p + k/q)/2, v = (k/q - j/p)/2, p = N + 1/2, q = N - 1/2\n")
                fh.write("# u and v are the generator parameters (fibre phase tau = pi/2 + 2*pi*N*u, i.e. a1 = pi/2 + 2*pi*p*u), not arclength fractions\n")
                fh.write("# sign = (-1)^k: A displaced to eta + sign*delta_eta and B to eta - sign*delta_eta at this crossing; x,y,z = crossing point on the undisplaced torus (unit ball)\n")
                fh.write("u,v,j,k,sign,x,y,z\n")
                for i in np.argsort(uc):
                    fh.write(f"{uc[i]:.9f},{vc[i]:.9f},{jj[i]},{kk[i]},{sg[i]},{X[i,0]:.6f},{X[i,1]:.6f},{X[i,2]:.6f}\n")
            num["files_extra"] = [f"mirror_pair_v1_{mode}.crossings.csv"]
        f2.mkdir(parents=True, exist_ok=True)
        key = f"mirror_pair_v1_{mode}"
        do_tube = tubes and mode == "woven"
        num["files"] = export_strands(f2, key, [("A", Af if do_tube else A, False, COL_A), ("B", Bf if do_tube else B, False, COL_B)], 0.0045, do_tube, 6,
                                      f"Hopf torus + z-mirror image, session-5 scheme ({mode}): eta {eta}, N {N}, revolutions 1, delta_eta {de}, tau0 pi/2, mirror rotation 0; "
                                      f"ccsim.geometry.hopf_mirror_pair({eta}, {N}, 1.0, '{mode}', {de}, tau0=np.pi/2, mirror_rotation=0.0)")
        if do_tube:
            write_obj_polylines(f2 / f"{key}.centreline.obj", [("A", A), ("B", B)], f"mirror pair v1 {mode} (canonical 3601-point polylines)")
            write_csv(f2 / f"{key}.A.csv", A, "mirror pair v1 strand A (canonical 3601 points)")
            write_csv(f2 / f"{key}.B.csv", B, "mirror pair v1 strand B = z-mirror of the undisplaced A with its own weave profile (canonical 3601 points)")
        manifest["geometries"][key] = num
        print(f"v1 {mode}: crossings={num.get('crossings')} minAB={dAB:.4f} altA={num.get('alternation_along_A')} altB={num.get('alternation_along_B')} |B-MA|={num['B_equals_zmirror_of_A_maxdev']:.4f}")

    # ======================================================================
    # 3. Meshed torus II — closed (N+1, N−1) knots, B = −A (rotated mirror = point inversion), zone weave, chiral option
    # ======================================================================
    f3 = OUT / "3_meshed_torus_II"
    cases = [("achiral_de015", 0.60, 0.0, 1, 0.015, 0.00448, True, "designs/woven_pair (η 0.60 variant): the plain mirror pair, version 2"),
             ("chiral_e070_n4_de015", 0.60, 0.70, 4, 0.015, 0.00448, True, "designs/woven_m24_de015.json — the design with nested flux surfaces (score 17.2 / 14.6 held-out)"),
             ("chiral_e070_n4_de010", 0.60, 0.70, 4, 0.010, 0.00303, False, "designs/woven_m24_de010.json — fully woven regime"),
             ("achiral_de010", 0.60, 0.0, 1, 0.010, 0.00303, False, "plain pair in the fully woven regime")]
    for key, eta, eps, n, de, a_wire, do_tube, note in cases:
        (A, B), info = G.hopf_helical_pair(eta, N, eps=eps, periods=n, delta_eta=de, return_crossings=True)
        sc = info["strand_scale"]
        eA, a1A, a2A = G.hopf_coordinates(A, sc)
        eB, a1B, a2B = G.hopf_coordinates(B, sc)
        dAB = cKDTree(A).query(B)[0].min()
        inv = float(np.abs(B + A).max())
        num = {"eta": eta, "circuits": N, "eps": eps, "periods": n, "delta_eta": de, "helical_mode": "toroidal", "mirror_rotation": "pi",
               "p": N + 1, "q": N - 1, "strand_scale": sc, "torus": frame(eta), "note": note,
               "A": strand_numbers(A, True), "B": strand_numbers(B, True), "min_distance_A_B": float(dAB),
               "B_equals_minus_A_maxdev": inv, "B_equals_minus_A": inv < 1e-9,
               "eta_range_A": [float(eA.min()), float(eA.max())], "eta_range_B": [float(eB.min()), float(eB.max())],
               "crossings": info["crossings"], "two_p_q": 2 * (N + 1) * (N - 1), "alternation_along_B": info["alternation_B"],
               "weave_regime": info["weave"], "zone_fraction": info["zone_fraction"], "weave_defects": info["weave_defects"],
               "components": info["components"], "min_distance_fine": info["min_distance_fine"], "closure_check": info["closure_check"],
               "conductor_radius_design": a_wire, "conductor_radius_used_for_tube": a_wire}
        # turn-density modulation check from the recovered coordinates: local da2/da1 along A vs 1 + eps cos(2 a2 - n a1)
        f3.mkdir(parents=True, exist_ok=True)
        uc, vc, sg = np.array(info["crossing_u"]), np.array(info["crossing_v"]), np.array(info["crossing_sign"])
        with open(f3 / f"meshed_torus_II_{key}.crossings.csv", "w", encoding="utf-8") as fh:
            fh.write("# crossings of the closed strand A(u) with B(v): u = a1/(2*pi*p) along A and v = a1/(2*pi*p) along the un-inverted B strand (p = N + 1; toroidal-angle fractions, not arclength);\n")
            fh.write("# sign = the over/under the zone weave used at this crossing (A displaced to eta + sign*delta_eta, B to eta - sign*delta_eta)\n")
            fh.write("u,v,sign\n")
            for i in np.argsort(uc):
                fh.write(f"{uc[i]:.9f},{vc[i]:.9f},{sg[i]}\n")
        num["files_extra"] = [f"meshed_torus_II_{key}.crossings.csv"]
        strands = [("A", A, True, COL_A), ("B", B, True, COL_B)]
        if do_tube:
            (Af, Bf), _ = G.hopf_helical_pair(eta, N, eps=eps, periods=n, delta_eta=de, n_per_circuit=200)
            strands_t = [("A", Af, True, COL_A), ("B", Bf, True, COL_B)]
        num["files"] = export_strands(f3, f"meshed_torus_II_{key}", strands_t if do_tube else strands, a_wire, do_tube, 6,
                                      f"Meshed torus II ({key}): eta {eta}, N {N}, eps {eps}, periods {n}, delta_eta {de}, closed (N+1, N-1) knots, B = R_z(pi)*M*A (= -A when eps = 0); "
                                      f"ccsim.geometry.hopf_helical_pair({eta}, {N}, eps={eps}, periods={n}, delta_eta={de})")
        if do_tube:
            write_obj_polylines(f3 / f"meshed_torus_II_{key}.centreline.obj", [("A", A), ("B", B)], f"meshed torus II {key} (canonical 3601-point closed polylines)")
            write_csv(f3 / f"meshed_torus_II_{key}.A.csv", A, f"meshed torus II {key} strand A (canonical 3601 points, closed: last = first)")
            write_csv(f3 / f"meshed_torus_II_{key}.B.csv", B, f"meshed torus II {key} strand B (canonical 3601 points, closed: last = first)")
        manifest["geometries"][f"meshed_torus_II_{key}"] = num
        print(f"v2 {key}: crossings={info['crossings']} weave={info['weave']} zones={info['zone_fraction']:.2f} minAB={dAB:.4f} |B+A|={inv:.2e} eta A {eA.min():.4f}..{eA.max():.4f}")

    manifest["seconds"] = time.time() - t_start
    json.dump(manifest, open(OUT / "MANIFEST.json", "w"), indent=1, default=float)
    print("manifest written;", time.time() - t_start, "s")


if __name__ == "__main__":
    main()
