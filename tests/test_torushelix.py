"""Geometry tests for the `torushelix` family (classical l=2 stellarator windings) — no numba needed.

    python3 -m pytest tests/test_torushelix.py -q      or      python3 tests/test_torushelix.py
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
from ccsim.design import FAMILIES, _parts, build, validate  # noqa: E402


def test_family_registered():
    assert "torushelix" in FAMILIES
    assert set(FAMILIES["torushelix"]["params"]) == {"R0", "r", "l", "periods", "phase_deg"}


def test_conductor_count_and_closure():
    for periods, expected in ((4, 4), (3, 2), (5, 2), (6, 4), (-4, 4)):
        parts = G.torus_helical_windings(0.6, 0.4, 2, periods)
        assert len(parts) == expected, (periods, len(parts))
        for pts, _ in parts:
            assert np.linalg.norm(pts[0] - pts[-1]) == 0.0          # exactly closed
            assert np.linalg.norm(pts, axis=1).max() <= 1.0 + 1e-12  # unit ball for R0 + r = 1


def test_signs_alternate_and_conductors_are_distinct():
    parts = G.torus_helical_windings(0.6, 0.4, 2, 4)
    assert [s for _, s in parts] == [1, -1, 1, -1]
    for i in range(len(parts)):
        for j in range(i + 1, len(parts)):
            d = np.linalg.norm(parts[i][0][:, None, :] - parts[j][0][None, ::7, :], axis=2).min()
            assert d > 0.1, (i, j, d)                                   # separate conductors, never coincident
    odd = G.torus_helical_windings(0.6, 0.4, 2, 3)
    assert [s for _, s in odd] == [1, -1]
    d = np.linalg.norm(odd[0][0][:, None, :] - odd[1][0][None, ::7, :], axis=2).min()
    assert d > 0.1


def test_on_torus_surface_and_handedness():
    R0, r = 0.6, 0.4
    (pts, _), = G.torus_helical_windings(R0, r, 2, 4)[:1]
    rho = np.hypot(pts[:, 0], pts[:, 1])
    assert np.allclose(np.hypot(rho - R0, pts[:, 2]), r)                # lies on the torus
    theta = np.unwrap(np.arctan2(pts[:, 2], rho - R0))
    phi = np.unwrap(np.arctan2(pts[:, 1], pts[:, 0]))
    slope = (theta[-1] - theta[0]) / (phi[-1] - phi[0])
    assert abs(slope - 2.0) < 1e-9                                        # dθ/dφ = n/l = 2
    (pts2, _), = G.torus_helical_windings(R0, r, 2, -4)[:1]
    theta2 = np.unwrap(np.arctan2(pts2[:, 2], np.hypot(pts2[:, 0], pts2[:, 1]) - R0))
    phi2 = np.unwrap(np.arctan2(pts2[:, 1], pts2[:, 0]))
    assert abs((theta2[-1] - theta2[0]) / (phi2[-1] - phi2[0]) + 2.0) < 1e-9   # opposite handedness


def test_design_parts_and_build():
    parts = _parts("torushelix", {})
    assert [(s, c) for _, s, c in parts] == [(1, True), (-1, True), (1, True), (-1, True)]
    d = {"name": "helix alone", "wire_radius": 0.005, "components": [{"family": "torushelix", "params": {"R0": 0.55, "r": 0.35, "l": 2, "periods": 4}, "current": 1.0}]}
    w = build(d)
    assert len(w.windings) == 4
    assert all(wd.active_count == len(wd.points) for wd in w.windings)     # closed: no remote return
    v = validate(d)
    assert v["valid"], v["problems"]


def test_composite_with_meshed_torus_validates():
    base = json.load(open(ROOT / "designs" / "woven_m24_de015.json"))
    d = {"name": "composite", "wire_radius": 0.0038,
         "components": [dict(base["components"][0], scale=0.85),
                        {"family": "torushelix", "params": {"R0": 0.543, "r": 0.365, "l": 2, "periods": 4}, "current": 0.6}]}
    v = validate(d)
    assert v["valid"], v["problems"]
    assert v["circuits"] == 6 and v["clearance"] > 2 * 0.0038


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
