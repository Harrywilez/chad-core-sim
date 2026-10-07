"""Declarative designs: components + transforms → a MultiWinding.

A design is a JSON-able dict:

    {"name": "two crossed precessing cores",
     "wire_radius": 0.005,                       # ball units (optional, default 0.005)
     "components": [
        {"family": "precess", "params": {"inward": 0.22, "rotation_deg": 24, "slip_deg": 3},
         "scale": 0.45,                          # uniform scale of the unit-ball generator
         "reflect": null,                        # null | "x" | "y" | "z"  (mirror plane normal)
         "rotate": {"axis": [1, 0, 0], "deg": 90},
         "translate": [0.5, 0, 0],
         "current": 1.0},                        # relative current (sign = direction)
        ...]}

Every family is one of the ccsim generators; a component may expand to several
circuits (mirror pairs, yin-yang, picket fence, links), each closed either on itself
(closed knots) or with a remote return.  Transform order: scale → reflect → rotate →
translate.  `FAMILIES` documents the parameter names and defaults an agent may use.
"""

from __future__ import annotations

import hashlib
import json
import math
from typing import Dict, List, Optional

import numpy as np

from . import geometry as G

FAMILIES: Dict[str, dict] = {
    "circle": {"params": {"radius": 0.7, "axis": "z"}, "doc": "single circular loop"},
    "solenoid": {"params": {"radius": 0.45, "length": 1.1, "turns": 12}, "doc": "straight solenoid along z"},
    "precess": {"params": {"inward": 0.22, "rotation_deg": 24.0, "slip_deg": 3.0}, "doc": "Codex precessing inward loops (mirror-type)"},
    "phase": {"params": {"inward": 0.32, "slip_deg": 3.0}, "doc": "Codex phase-slip toroidal spiral"},
    "codexhopf": {"params": {"inward": 0.18, "rotation_deg": 20.0}, "doc": "Codex Hopf-coordinate drift trial"},
    "hopfdrift": {"params": {"variant": "s64_180", "circuits": 10}, "doc": "Hopf drift winding; variants s64_150/180/210, d48_150/180/210"},
    "hopfcont": {"params": {"eta0": 0.82, "eta1": 0.64, "sweep_deg": 720, "circuits": 48}, "doc": "Hopf drift continued toward a torus"},
    "hopftorus": {"params": {"eta": 0.70, "circuits": 24, "revolutions": 1.0}, "doc": "single Hopf-torus strand ((N+r/2, N−r/2) torus curve)"},
    "recursive": {"params": {"base_radius": 0.72, "levels": [[0.24, 12]], "alternate": False}, "doc": "recursive ∫∫dl / ∫∫∫dl winding; levels = [[rel_radius, turns], ...]"},
    "baseball": {"params": {"radius": 0.72, "amplitude_deg": 40.0, "turns": 6, "radial_pitch": 0.035}, "doc": "baseball-seam (min-B) coil"},
    "sphere": {"params": {"radius": 0.76, "turns": 24}, "doc": "sphere winding, uniform interior field"},
    "hopfmirror": {"params": {"eta": 0.70, "circuits": 24, "revolutions": 2.0, "mode": "woven", "delta_eta": 0.03, "sense": -1}, "doc": "Hopf torus + exact rotated mirror image (2 circuits); sense −1 = opposed currents (toroidal solenoid), +1 = same"},
    "hopfmesh": {"params": {"eta": 0.60, "circuits": 24, "eps": 0.70, "periods": 4, "delta_eta": 0.03, "helical_mode": "toroidal", "sense": -1}, "doc": "chiral meshed torus (closed knots, helical density modulation) — the stellarator weave"},
    "yinyang": {"params": {"radius_outer": 0.76, "radius_inner": 0.60, "amplitude_deg": 40.0, "turns": 4, "sense": 1}, "doc": "two interleaved baseball coils (2 circuits)"},
    "picket": {"params": {"rings": 5, "radius": 0.72, "alternating": True}, "doc": "ring-cusp stack (one circuit per ring)"},
    "link": {"params": {"radius": 0.62, "sense": 1}, "doc": "two linked rings (Hopf link)"},
    "cone": {"params": {"radius_base": 0.32, "radius_tip": 0.06, "height": 0.5, "turns": 8}, "doc": "'tornado' coil: conical helix along z, base at −height/2, tip (small radius, strongest field) at +height/2; reflect z to point the tip down"},
    "torushelix": {"params": {"R0": 0.60, "r": 0.40, "l": 2, "periods": 4, "phase_deg": 0.0}, "doc": "classical stellarator helical windings on the torus (R0, r): 2l closed conductors θ = (n/l)φ + kπ/l with alternating currents (2l circuits); negative periods = opposite handedness; combine with a toroidal-field component"},
}


def _parts(family: str, p: dict) -> List[tuple]:
    """Unit-ball point sets for a family: list of (points, sign, closed_on_itself)."""
    d = dict(FAMILIES[family]["params"]); d.update(p or {})
    if family == "circle":
        return [(G.circular_loop(d["radius"], 720, axis=d["axis"]), 1, True)]
    if family == "solenoid":
        return [(G.solenoid(d["radius"], d["length"], int(d["turns"]), n_per_turn=60), 1, False)]
    if family == "precess":
        return [(G.codex_toroidal("precess", d["inward"], d["rotation_deg"], d["slip_deg"]), 1, False)]
    if family == "phase":
        return [(G.codex_toroidal("phase", d["inward"], 3.0, d["slip_deg"]), 1, False)]
    if family == "codexhopf":
        return [(G.codex_hopf(d["inward"], d["rotation_deg"]), 1, False)]
    if family == "hopfdrift":
        return [(G.hopf_drift(d["variant"], int(d["circuits"])), 1, False)]
    if family == "hopfcont":
        return [(G.hopf_continued(d["eta0"], d["eta1"], d["sweep_deg"], int(d["circuits"])), 1, False)]
    if family == "hopftorus":
        pts = G.hopf_torus(d["eta"], int(d["circuits"]), d["revolutions"])
        closed = abs(d["revolutions"] - round(d["revolutions"])) < 1e-9 and int(round(d["revolutions"])) % 2 == 0
        return [(pts, 1, closed)]
    if family == "recursive":
        lv = [(float(a), int(b)) for a, b in d["levels"]]
        return [(G.recursive_winding("circle", d["base_radius"], levels=lv, n_per_turn=24 if len(lv) > 1 else 48, alternate_handedness=bool(d["alternate"])), 1, False)]
    if family == "baseball":
        return [(G.baseball_seam(d["radius"], d["amplitude_deg"], int(d["turns"]), d["radial_pitch"]), 1, False)]
    if family == "sphere":
        return [(G.sphere_winding(d["radius"], int(d["turns"])), 1, False)]
    if family == "hopfmirror":
        if abs(d["revolutions"] - 2.0) < 1e-9 and d["mode"] != "nested":
            # closed (N+1, N−1) knots: the unmodulated meshed torus — same generator (and zone weave) as hopfmesh with ε = 0
            (A, B), _ = G.hopf_helical_pair(d["eta"], int(d["circuits"]), eps=0.0, periods=1, delta_eta=d["delta_eta"])
            return [(A, 1, True), (B, int(d["sense"]), True)]
        A, B = G.hopf_mirror_pair(d["eta"], int(d["circuits"]), d["revolutions"], d["mode"], d["delta_eta"])
        closed = np.linalg.norm(A[0] - A[-1]) < 1e-9
        return [(A, 1, closed), (B, int(d["sense"]), closed)]
    if family == "hopfmesh":
        (A, B), _ = G.hopf_helical_pair(d["eta"], int(d["circuits"]), eps=d["eps"], periods=int(d["periods"]), delta_eta=d["delta_eta"], helical_mode=d["helical_mode"])
        return [(A, 1, True), (B, int(d["sense"]), True)]
    if family == "yinyang":
        o, i = G.yin_yang(d["radius_outer"], d["radius_inner"], d["amplitude_deg"], int(d["turns"]))
        return [(o, 1, False), (i, int(d["sense"]), False)]
    if family == "picket":
        return [(pts, sgn if d["alternating"] else 1, True) for pts, sgn in G.picket_fence(int(d["rings"]), d["radius"])]
    if family == "link":
        a, b = G.hopf_link(d["radius"])
        return [(a, 1, True), (b, int(d["sense"]), True)]
    if family == "cone":
        return [(G.cone_helix(d["radius_base"], d["radius_tip"], d["height"], int(d["turns"])), 1, False)]
    if family == "torushelix":
        return [(pts, sgn, True) for pts, sgn in G.torus_helical_windings(d["R0"], d["r"], int(d["l"]), int(d["periods"]), d["phase_deg"])]
    raise KeyError(f"unknown family {family!r}; known: {sorted(FAMILIES)}")


def transform(pts: np.ndarray, comp: dict) -> np.ndarray:
    out = np.asarray(pts, dtype=float) * float(comp.get("scale", 1.0))
    refl = comp.get("reflect")
    if refl:
        m = np.ones(3); m["xyz".index(refl)] = -1.0
        out = out * m
    rot = comp.get("rotate")
    if rot:
        Rm = G.rotation_matrix(rot.get("axis", [0, 0, 1]), math.radians(rot.get("deg", 0.0)))
        out = out @ Rm.T
    tr = comp.get("translate")
    if tr:
        out = out + np.asarray(tr, dtype=float)
    return out


def build(design: dict, return_radius_factor: float = 3.8) -> "G.MultiWinding":
    """Design dict → MultiWinding (unit ball, unit current per component × `current`)."""
    a = float(design.get("wire_radius", 0.005))
    ws = []
    for ci, comp in enumerate(design["components"]):
        fam = comp["family"]
        if float(comp.get("current", 1.0)) == 0.0:      # a zero-current component is disabled (no conductor, no field)
            continue
        for k, (pts, sgn, closed) in enumerate(_parts(fam, comp.get("params", {}))):
            P = transform(pts, comp)
            cur = float(comp.get("current", 1.0)) * sgn
            if closed and np.linalg.norm(P[0] - P[-1]) < 1e-9:
                closed_pts, n_active = P, len(P)
            else:
                closed_pts, n_active = G.close_with_return(P, return_radius_factor)
            ws.append(G.Winding(f"{design.get('name', 'design')} [{ci}:{fam} #{k}]", closed_pts, cur, a, closed=True, active_count=n_active,
                                meta={"component": ci, "family": fam, "part": k, "sign": cur}))
    if not ws:
        raise ValueError("design has no component with a non-zero current")
    return G.MultiWinding(design.get("name", "design"), ws, {"family": "design", "design": design})


def design_hash(design: dict) -> str:
    return hashlib.sha1(json.dumps(design, sort_keys=True, default=float).encode()).hexdigest()[:12]


def validate(design: dict, min_clearance: Optional[float] = None) -> dict:
    """Geometric sanity: inside the unit ball, conductors not touching (clearance ≥ 2 × wire
    radius unless `min_clearance` is given), size of the job."""
    if min_clearance is None:
        min_clearance = 2.0 * float(design.get("wire_radius", 0.005))
    w = build(design)
    act = w.active_points
    rmax = float(np.linalg.norm(act, axis=1).max())
    clr = float(w.clearance_m())
    npts = int(sum(len(x.points) for x in w.windings))
    problems = []
    if rmax > 1.0 + 1e-6:
        problems.append(f"conductor reaches r = {rmax:.3f} > 1 (outside the ball)")
    if clr < min_clearance:
        problems.append(f"conductor clearance {clr:.4f} < 2 × wire radius {min_clearance:.4f} (conductors touch; thin the wire or open the geometry)")
    if npts > 60000:
        problems.append(f"{npts} conductor points — too expensive; reduce turns/circuits")
    return {"valid": not problems, "problems": problems, "r_max": rmax, "clearance": clr, "conductor_points": npts, "circuits": len(w.windings), "hash": design_hash(design)}
