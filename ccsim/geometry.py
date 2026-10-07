"""Closed-wire winding geometries.

All generators return an (N, 3) polyline in **metres** describing a single
conductor.  A `Winding` also carries the conductor radius and the current, so
the field modules can treat it as a complete magnetostatic source.

Three families are provided:

1. **Recursive windings** — the Chad Core Prime concept as sketched:
   a wire at 1× is a line, at 10× it is a coil, at 100× it is a coil of coils,
   "and so on".  `recursive_winding` wraps a helix of radius r_k with n_k turns
   around the level-(k−1) curve using a rotation-minimising (Bishop) frame, so
   any depth is possible and straight parent segments are handled.
2. **Codex realizations** (ported from the 2026-08-31 closed-wire trials and
   Hopf factorial so the numbers are comparable): phase-slip toroidal spiral,
   precessing inward loops, Hopf-coordinate drift, and the six Hopf-drift
   (eta1, sweep) variants.  They are generated in the unit ball and scaled.
3. **Reference coils** for verification: circular loop, finite solenoid,
   Helmholtz pair, straight wire.

Every open winding can be closed with `remote_return`, which adds a smooth
return conductor well outside the vessel; a closed circuit is what makes
∇·J = 0 hold and Ampère's law checkable.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from scipy.spatial import cKDTree
from typing import Optional, Sequence, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# Basic polyline utilities
# ---------------------------------------------------------------------------


def unit(v: np.ndarray) -> np.ndarray:
    v = np.asarray(v, dtype=float)
    n = np.linalg.norm(v, axis=-1, keepdims=True)
    return v / np.maximum(n, 1e-300)


def resample_arclength(points: np.ndarray, count: int) -> np.ndarray:
    """Resample a polyline to `count` points equally spaced in arclength."""

    p = np.asarray(points, dtype=float)
    ds = np.linalg.norm(np.diff(p, axis=0), axis=1)
    keep = np.r_[True, ds > 1e-14]
    p = p[keep]
    s = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(p, axis=0), axis=1))]
    q = np.linspace(0.0, s[-1], count)
    return np.column_stack([np.interp(q, s, p[:, j]) for j in range(3)])


def polyline_length(points: np.ndarray) -> float:
    return float(np.sum(np.linalg.norm(np.diff(np.asarray(points), axis=0), axis=1)))


def smoothstep5(u: np.ndarray) -> np.ndarray:
    return u**3 * (10.0 + u * (-15.0 + 6.0 * u))


def bishop_frame(points: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Rotation-minimising frame (T, N1, N2) along a polyline (double reflection)."""

    p = np.asarray(points, dtype=float)
    n = len(p)
    t = np.zeros_like(p)
    t[:-1] = unit(p[1:] - p[:-1])
    t[-1] = t[-2]
    # smooth tangents by averaging adjacent segment directions
    t[1:-1] = unit(t[:-2] + t[1:-1])
    n1 = np.zeros_like(p)
    # initial normal: any vector perpendicular to t[0]
    trial = np.array([0.0, 0.0, 1.0]) if abs(t[0][2]) < 0.9 else np.array([1.0, 0.0, 0.0])
    n1[0] = unit(trial - np.dot(trial, t[0]) * t[0])
    for i in range(n - 1):
        v1 = p[i + 1] - p[i]
        c1 = float(np.dot(v1, v1))
        if c1 < 1e-30:
            n1[i + 1] = n1[i]
            continue
        rL = n1[i] - (2.0 / c1) * np.dot(v1, n1[i]) * v1
        tL = t[i] - (2.0 / c1) * np.dot(v1, t[i]) * v1
        v2 = t[i + 1] - tL
        c2 = float(np.dot(v2, v2))
        n1[i + 1] = rL if c2 < 1e-30 else rL - (2.0 / c2) * np.dot(v2, rL) * v2
        n1[i + 1] = unit(n1[i + 1] - np.dot(n1[i + 1], t[i + 1]) * t[i + 1])
    n2 = np.cross(t, n1)
    return t, n1, n2


def min_nonlocal_distance(points: np.ndarray, skip: int = 10) -> Tuple[float, Tuple[int, int]]:
    """Closest approach between non-adjacent parts of a polyline (conductor clearance)."""

    p = np.asarray(points, dtype=float)
    n = len(p)
    best = float("inf")
    pair = (0, skip + 1)
    block = 128
    for i0 in range(0, n, block):
        i1 = min(n, i0 + block)
        a = p[i0:i1]
        d2 = np.sum((a[:, None, :] - p[None, :, :]) ** 2, axis=2)
        rows = np.arange(i0, i1)[:, None]
        cols = np.arange(n)[None, :]
        sep = np.abs(rows - cols)
        d2[(sep <= skip) | (sep >= n - skip)] = np.inf
        local = int(np.argmin(d2))
        val = float(d2.flat[local])
        if val < best * best:
            rr, cc = np.unravel_index(local, d2.shape)
            best = math.sqrt(val)
            pair = (i0 + int(rr), int(cc))
    return best, pair


def nearest_distance(points: np.ndarray, samples: np.ndarray) -> np.ndarray:
    """Distance from each point to the nearest sample point (brute force, chunked)."""

    pts = np.atleast_2d(points)
    result = np.full(len(pts), np.inf)
    for i in range(0, len(samples), 512):
        d2 = np.sum((pts[:, None, :] - samples[None, i : i + 512, :]) ** 2, axis=2)
        result = np.minimum(result, np.sqrt(np.min(d2, axis=1)))
    return result


# ---------------------------------------------------------------------------
# Winding container
# ---------------------------------------------------------------------------


@dataclass
class Winding:
    """A single closed or open conductor path with its current and radius."""

    name: str
    points: np.ndarray          # (N, 3) metres
    current_A: float            # amperes along the path direction
    wire_radius_m: float
    closed: bool = False
    active_count: Optional[int] = None   # number of points belonging to the "active" (in-vessel) part
    meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.points = np.asarray(self.points, dtype=float)
        if self.points.ndim != 2 or self.points.shape[1] != 3 or len(self.points) < 2:
            raise ValueError("points must be an (N>=2, 3) array")
        if self.wire_radius_m <= 0:
            raise ValueError("wire radius must be positive")
        if self.active_count is None:
            self.active_count = len(self.points)

    @property
    def active_points(self) -> np.ndarray:
        return self.points[: self.active_count]

    @property
    def length_m(self) -> float:
        return polyline_length(self.points)

    @property
    def active_length_m(self) -> float:
        return polyline_length(self.active_points)

    def segments(self) -> Tuple[np.ndarray, np.ndarray]:
        """Segment start points and vectors (both (N-1, 3))."""

        return self.points[:-1], np.diff(self.points, axis=0)

    def clearance_m(self, skip: int = 12) -> float:
        return min_nonlocal_distance(self.active_points, skip=skip)[0]

    def scaled(self, factor: float, name: Optional[str] = None) -> "Winding":
        return Winding(name or self.name, self.points * factor, self.current_A, self.wire_radius_m * factor,
                       self.closed, self.active_count, dict(self.meta))

    def with_current(self, current_A: float) -> "Winding":
        return Winding(self.name, self.points, current_A, self.wire_radius_m, self.closed, self.active_count, dict(self.meta))

    def with_wire_radius(self, wire_radius_m: float) -> "Winding":
        return Winding(self.name, self.points, self.current_A, wire_radius_m, self.closed, self.active_count, dict(self.meta))

    def resistance_ohm(self, resistivity_ohm_m: float = 1.68e-8) -> float:
        """DC resistance of the whole path (copper by default)."""

        return resistivity_ohm_m * self.length_m / (math.pi * self.wire_radius_m**2)

    def joule_power_W(self, resistivity_ohm_m: float = 1.68e-8) -> float:
        return self.current_A**2 * self.resistance_ohm(resistivity_ohm_m)


# ---------------------------------------------------------------------------
# Reference coils
# ---------------------------------------------------------------------------


def circular_loop(radius_m: float, n: int = 720, center=(0.0, 0.0, 0.0), axis="z") -> np.ndarray:
    th = np.linspace(0.0, 2.0 * np.pi, n + 1)
    c = np.asarray(center, dtype=float)
    if axis == "z":
        return c + np.column_stack((radius_m * np.cos(th), radius_m * np.sin(th), np.zeros_like(th)))
    if axis == "x":
        return c + np.column_stack((np.zeros_like(th), radius_m * np.cos(th), radius_m * np.sin(th)))
    return c + np.column_stack((radius_m * np.sin(th), np.zeros_like(th), radius_m * np.cos(th)))


def solenoid(radius_m: float, length_m: float, turns: int, n_per_turn: int = 90) -> np.ndarray:
    n = turns * n_per_turn + 1
    th = np.linspace(0.0, 2.0 * np.pi * turns, n)
    z = np.linspace(-0.5 * length_m, 0.5 * length_m, n)
    return np.column_stack((radius_m * np.cos(th), radius_m * np.sin(th), z))


def straight_wire(length_m: float, n: int = 2) -> np.ndarray:
    z = np.linspace(-0.5 * length_m, 0.5 * length_m, n)
    return np.column_stack((np.zeros_like(z), np.zeros_like(z), z))


def helmholtz_pair(radius_m: float, n: int = 720) -> np.ndarray:
    """Two coaxial loops separated by R, joined by a short axial jumper (single conductor)."""

    a = circular_loop(radius_m, n, center=(0, 0, -0.5 * radius_m))
    b = circular_loop(radius_m, n, center=(0, 0, 0.5 * radius_m))
    # go up to the second loop, around it, and back down the same line: the two
    # jumper segments are antiparallel and cancel exactly, leaving two pure loops.
    return np.vstack((a, b, a[:1]))


# ---------------------------------------------------------------------------
# Recursive windings: the ∫dl, ∫∫dl, ∫∫∫dl sketch
# ---------------------------------------------------------------------------


def helix_around(parent: np.ndarray, radius_m: float, turns: float, n_per_turn: int = 48,
                 phase: float = 0.0, handedness: float = 1.0) -> np.ndarray:
    """Wrap a helix of `turns` turns and radius `radius_m` around a parent polyline.

    The parent is first resampled in arclength so the winding pitch is uniform;
    a Bishop frame avoids the Frenet singularity on straight or planar parents.
    """

    n = max(int(turns * n_per_turn) + 1, 8)
    base = resample_arclength(parent, n)
    t, n1, n2 = bishop_frame(base)
    ang = phase + handedness * 2.0 * np.pi * turns * np.linspace(0.0, 1.0, n)
    return base + radius_m * (np.cos(ang)[:, None] * n1 + np.sin(ang)[:, None] * n2)


def recursive_winding(base: str = "circle", base_radius_m: float = 0.10, base_length_m: float = 0.30,
                      levels: Sequence[Tuple[float, float]] = ((0.3, 12),),
                      n_per_turn: int = 48, alternate_handedness: bool = False) -> np.ndarray:
    """Build the sketch's nested winding.

    `levels` is a sequence of (relative radius, turns) pairs: level-1 has
    radius levels[0][0] * base_radius and levels[0][1] turns around the base
    curve; level-2 has radius levels[1][0] * (level-1 radius) and levels[1][1]
    turns *per level-1 turn*, and so on.  A wire at 1× is the base ('line' or
    'circle'); at 10× the first helix; at 100× the second.
    """

    if base == "circle":
        turns0 = float(levels[0][1]) if levels else 12.0
        th = np.linspace(0.0, 2.0 * np.pi * (1.0 - 0.5 / turns0), max(360, 24 * int(turns0)))
        curve = np.column_stack((base_radius_m * np.cos(th), base_radius_m * np.sin(th), np.zeros_like(th)))
    elif base == "line":
        curve = straight_wire(base_length_m, n=64)
    else:
        raise ValueError("base must be 'circle' or 'line'")
    radius = base_radius_m
    turns_total = 1.0
    for k, (rel_radius, turns_per_parent_turn) in enumerate(levels):
        radius = rel_radius * radius
        turns_total *= turns_per_parent_turn
        hand = -1.0 if (alternate_handedness and k % 2 == 1) else 1.0
        curve = helix_around(curve, radius, turns_total, n_per_turn=n_per_turn, handedness=hand)
    return curve


# ---------------------------------------------------------------------------
# Codex realizations (ported; generated in the unit ball, then scaled)
# ---------------------------------------------------------------------------

_N_CIRCUITS = 10
_K_MINOR = 11
_ACTIVE_SEGMENTS = 1800


def _rotate_z(v, angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.column_stack((c * v[:, 0] - s * v[:, 1], s * v[:, 0] + c * v[:, 1], v[:, 2]))


def _rotate_y(v, angle):
    c, s = math.cos(angle), math.sin(angle)
    return np.column_stack((c * v[:, 0] + s * v[:, 2], v[:, 1], -s * v[:, 0] + c * v[:, 2]))


def codex_toroidal(kind: str = "phase", inward: float = 0.32, rotation_deg: float = 3.0,
                   minor_slip_deg: float = 3.0, n_hi: int = 8193) -> np.ndarray:
    """Phase-slip ('phase') or precessing ('precess') inward toroidal spiral, unit ball."""

    u = np.linspace(0.0, 1.0, n_hi)
    h = smoothstep5(u)
    loops = _N_CIRCUITS * u
    theta = 2.0 * np.pi * loops
    r_outer = 4.45
    r_inner = 4.45 - inward * _N_CIRCUITS
    radius = r_outer + (r_inner - r_outer) * h
    envelope = 0.52 + (0.52 * 0.88 - 0.52) * h
    psi = _K_MINOR * theta + np.deg2rad(minor_slip_deg) * loops
    v = np.column_stack(
        ((radius + envelope * np.cos(psi)) * np.cos(theta),
         (radius + envelope * np.cos(psi)) * np.sin(theta),
         envelope * np.sin(psi))
    )
    if kind == "precess":
        chi = np.deg2rad(rotation_deg) * loops
        v = _rotate_z(_rotate_y(_rotate_z(v, -chi), math.radians(28.0)), chi)
    p = resample_arclength(v, _ACTIVE_SEGMENTS + 1)
    return p / np.linalg.norm(p, axis=1).max()


def codex_hopf(inward: float = 0.18, rotation_deg: float = 20.0, n_hi: int = 8193) -> np.ndarray:
    """Hopf-coordinate drift winding from the closed-wire trials, unit ball."""

    u = np.linspace(0.0, 1.0, n_hi)
    h = smoothstep5(u)
    loops = _N_CIRCUITS * u
    tau = 2.0 * np.pi * loops
    eta0 = 0.82
    eta1 = max(0.34, eta0 - (0.16 + 0.17 * inward))
    eta = eta0 + (eta1 - eta0) * h
    phi = 0.18 + np.deg2rad(rotation_deg) * loops
    a1 = tau + phi / 2.0
    a2 = tau - phi / 2.0
    q1 = np.cos(eta) * np.cos(a1)
    q2 = np.cos(eta) * np.sin(a1)
    q3 = np.sin(eta) * np.cos(a2)
    q4 = np.sin(eta) * np.sin(a2)
    den = 1.0 - q4
    if den.min() <= 1e-3:
        raise ValueError("Hopf stereographic projection approached its pole")
    p = np.column_stack((q1 / den, q2 / den, q3 / den))
    p = resample_arclength(p, _ACTIVE_SEGMENTS + 1)
    return p / np.linalg.norm(p, axis=1).max()


HOPF_DRIFT_VARIANTS = {
    "s64_150": {"eta1": 0.64, "sweep": 150.0, "label": "Shallow 150°"},
    "s64_180": {"eta1": 0.64, "sweep": 180.0, "label": "Shallow 180°"},
    "s64_210": {"eta1": 0.64, "sweep": 210.0, "label": "Shallow 210°"},
    "d48_150": {"eta1": 0.48, "sweep": 150.0, "label": "Deep 150°"},
    "d48_180": {"eta1": 0.48, "sweep": 180.0, "label": "Deep 180°"},
    "d48_210": {"eta1": 0.48, "sweep": 210.0, "label": "Deep 210°"},
}


def hopf_drift(variant: str = "s64_180", circuits: int = 10, n_per_circuit: int = 150) -> np.ndarray:
    """Hopf-drift winding from the factorial screen (unit ball, arclength-uniform eta sweep).

    Fibre phase τ = π/2 + 2π·circuits·u; the fibre (η, φ) drifts from (0.82, 0) to (η₁, sweep) at
    constant speed on the Hopf base sphere (metric dη² + ¼ sin²2η dφ²); scaled so the outer equator
    of the starting torus η₀ sits at radius 1 (scale (1 − sin η₀)/cos η₀)."""

    spec = HOPF_DRIFT_VARIANTS[variant]
    n_hi = max(4097, circuits * 512 + 1)
    u = np.linspace(0.0, 1.0, n_hi)
    eta0 = 0.82
    delta_eta = spec["eta1"] - eta0
    delta_phi = np.deg2rad(spec["sweep"])
    v_table = np.linspace(0.0, 1.0, 4097)
    eta_table = eta0 + delta_eta * v_table
    speed = np.sqrt(4.0 * delta_eta**2 + np.sin(2.0 * eta_table) ** 2 * delta_phi**2)
    cumulative = np.r_[0.0, np.cumsum(0.5 * (speed[:-1] + speed[1:]) * np.diff(v_table))]
    cumulative /= cumulative[-1]
    v = np.interp(u, cumulative, v_table)
    eta = eta0 + delta_eta * v
    phi = delta_phi * v
    tau = np.pi / 2.0 + 2.0 * np.pi * circuits * u
    a1 = tau + phi / 2.0
    a2 = tau - phi / 2.0
    q1 = np.cos(eta) * np.cos(a1)
    q2 = np.cos(eta) * np.sin(a1)
    q3 = np.sin(eta) * np.cos(a2)
    q4 = np.sin(eta) * np.sin(a2)
    den = 1.0 - q4
    if den.min() <= 1e-3:
        raise ValueError(f"{variant} approached the stereographic pole")
    scale = (1.0 - np.sin(eta0)) / np.cos(eta0)
    points = scale * np.column_stack((q1 / den, q2 / den, q3 / den))
    return resample_arclength(points, circuits * n_per_circuit + 1)


# ---------------------------------------------------------------------------
# Closing the circuit
# ---------------------------------------------------------------------------


def _bezier(p0, p1, p2, p3, n):
    t = np.linspace(0.0, 1.0, n)[:, None]
    return (1 - t) ** 3 * p0 + 3 * (1 - t) ** 2 * t * p1 + 3 * (1 - t) * t**2 * p2 + t**3 * p3


def remote_return(active: np.ndarray, radius: float, n_points: int = 300) -> np.ndarray:
    """Smooth return conductor from the end of `active` back to its start, routed
    outside a sphere of the given radius (metres).  Ported from the Codex trials
    so return-path effects stay comparable; use `return_sensitivity` in
    fields.py to quantify how much the return contributes inside the vessel."""

    start, end = active[0], active[-1]
    t0 = unit(active[1] - active[0])
    t1 = unit(active[-1] - active[-2])
    n0, n1 = unit(start), unit(end)
    q0, q1 = radius * n0, radius * n1
    axis = np.cross(n1, n0)
    if np.linalg.norm(axis) < 0.12:
        trial = np.array([0.0, 0.0, 1.0])
        if abs(np.dot(trial, n0)) > 0.85:
            trial = np.array([0.0, 1.0, 0.0])
        axis = np.cross(n0, trial)
    axis = unit(axis)
    if axis[2] < 0:
        axis = -axis
    apex = (radius * 1.24) * axis
    tangent = unit(q0 - q1)
    leg = 0.09 * radius
    outer = 0.22 * radius
    parts = [
        _bezier(end, end + leg * t1, q1 - outer * n1, q1, 80),
        _bezier(q1, q1 + outer * n1, apex - outer * tangent, apex, 95),
        _bezier(apex, apex + outer * tangent, q0 + outer * n0, q0, 95),
        _bezier(q0, q0 - outer * n0, start - leg * t0, start, 80),
    ]
    raw = np.vstack([parts[0], parts[1][1:], parts[2][1:], parts[3][1:]])
    return resample_arclength(raw, n_points + 1)


def close_with_return(active: np.ndarray, return_radius_m: float, n_return: int = 300) -> Tuple[np.ndarray, int]:
    """Return (closed polyline, active_count)."""

    ret = remote_return(active, return_radius_m, n_return)
    return np.vstack((active, ret[1:])), len(active)


# ---------------------------------------------------------------------------
# Catalogue of Chad-core configurations
# ---------------------------------------------------------------------------


def chadcore_catalogue(scale_m: float, wire_radius_m: float, current_A: float,
                       return_radius_factor: float = 3.8, include_hopf_drift: bool = True) -> list:
    """All Chad-core configurations as closed Windings at a given vessel scale.

    `scale_m` is the radius of the unit ball the Codex windings are generated in
    (their wall sits at 0.82 of it).  Recursive windings are sized so that
    their outermost extent matches the same ball.
    """

    items = []

    def add(name, active, meta):
        pts, n_active = close_with_return(active, return_radius_factor * scale_m)
        items.append(Winding(name, pts, current_A, wire_radius_m, closed=True, active_count=n_active, meta=meta))

    # --- recursive (∫dl → ∫∫dl → ∫∫∫dl) ---
    lvl1 = recursive_winding("circle", base_radius_m=0.72 * scale_m, levels=[(0.24, 12)])
    add("recursive-L2 (loop of loops, 12 turns)", lvl1, {"family": "recursive", "levels": 2})
    lvl2 = recursive_winding("circle", base_radius_m=0.72 * scale_m, levels=[(0.24, 12), (0.30, 6)], n_per_turn=24)
    add("recursive-L3 (coil of coils, 12×6)", lvl2, {"family": "recursive", "levels": 3})
    lvl2b = recursive_winding("circle", base_radius_m=0.72 * scale_m, levels=[(0.24, 12), (0.30, 6)], n_per_turn=24,
                              alternate_handedness=True)
    add("recursive-L3 alt-hand (12×6, opposite handedness)", lvl2b, {"family": "recursive", "levels": 3, "alternate": True})

    # --- Codex realizations (unit ball × scale) ---
    add("codex phase-slip toroidal spiral", scale_m * codex_toroidal("phase", 0.32, 3.0, 3.0),
        {"family": "codex", "kind": "phase"})
    add("codex precessing inward loops", scale_m * codex_toroidal("precess", 0.22, 24.0, 3.0),
        {"family": "codex", "kind": "precess"})
    add("codex Hopf-coordinate drift", scale_m * codex_hopf(0.18, 20.0), {"family": "codex", "kind": "hopf"})

    if include_hopf_drift:
        for key, spec in HOPF_DRIFT_VARIANTS.items():
            add(f"hopf-drift {spec['label']} (10 circuits)", scale_m * hopf_drift(key, 10),
                {"family": "hopf-drift", "variant": key, "circuits": 10})

    # --- reference: a plain solenoid of comparable size (control) ---
    sol = solenoid(0.45 * scale_m, 1.1 * scale_m, 12, n_per_turn=60)
    add("reference solenoid (12 turns)", sol, {"family": "reference"})
    return items


# ---------------------------------------------------------------------------
# Hopf continuation: keep the fibre drifting until the winding fills a torus
# ---------------------------------------------------------------------------


def _hopf_points(eta: np.ndarray, phi: np.ndarray, tau: np.ndarray) -> np.ndarray:
    """Stereographic image of the S³ curve (cosη e^{i a1}, sinη e^{i a2}), a1 = τ+φ/2, a2 = τ−φ/2."""

    a1 = tau + phi / 2.0
    a2 = tau - phi / 2.0
    q1 = np.cos(eta) * np.cos(a1)
    q2 = np.cos(eta) * np.sin(a1)
    q3 = np.sin(eta) * np.cos(a2)
    q4 = np.sin(eta) * np.sin(a2)
    den = 1.0 - q4
    if den.min() <= 1e-3:
        raise ValueError("stereographic pole reached; keep eta < ~0.9")
    return np.column_stack((q1 / den, q2 / den, q3 / den))


def hopf_coordinates(points: np.ndarray, scale: float = 1.0) -> tuple:
    """Inverse of `_hopf_points`: (η, a1, a2) of points given in the generator's *raw* frame divided by
    `scale` (pass the strand_scale reported by the pair generators for normalised output).  Inverse
    stereographic projection q = (2p, |p|² − 1)/(|p|² + 1), then η = atan2(|(q3, q4)|, |(q1, q2)|),
    a1 = atan2(q2, q1), a2 = atan2(q4, q3).  Toroidal angle φ_geo = a1; the fibre angle is
    τ = (a1 + a2)/2 and the base-sphere azimuth φ = a1 − a2."""

    p = np.asarray(points, dtype=float) / scale
    r2 = np.sum(p * p, axis=1)
    q = np.column_stack((2 * p / (r2 + 1)[:, None], (r2 - 1) / (r2 + 1)))
    eta = np.arctan2(np.hypot(q[:, 2], q[:, 3]), np.hypot(q[:, 0], q[:, 1]))
    return eta, np.arctan2(q[:, 1], q[:, 0]), np.arctan2(q[:, 3], q[:, 2])


def hopf_continued(eta0: float = 0.82, eta1: float = 0.64, sweep_deg: float = 720.0, circuits: int = 48,
                   n_per_circuit: int = 150) -> np.ndarray:
    """The Hopf-drift winding continued: same fibre drift law, but φ keeps sweeping
    (720° = twice round the torus) while η eases from eta0 to eta1, so the
    conductor covers the whole toroidal shell instead of a 150° wedge.  Unit ball."""

    n_hi = max(4097, circuits * 400 + 1)
    u = np.linspace(0.0, 1.0, n_hi)
    delta_eta = eta1 - eta0
    delta_phi = np.deg2rad(sweep_deg)
    v_table = np.linspace(0.0, 1.0, 4097)
    eta_table = eta0 + delta_eta * v_table
    speed = np.sqrt(4.0 * delta_eta**2 + np.sin(2.0 * eta_table) ** 2 * delta_phi**2)
    cumulative = np.r_[0.0, np.cumsum(0.5 * (speed[:-1] + speed[1:]) * np.diff(v_table))]
    cumulative /= cumulative[-1]
    v = np.interp(u, cumulative, v_table)
    eta = eta0 + delta_eta * v
    phi = delta_phi * v
    tau = np.pi / 2.0 + 2.0 * np.pi * circuits * u
    pts = _hopf_points(eta, phi, tau)
    pts *= (1.0 - np.sin(eta0)) / np.cos(eta0)
    pts = resample_arclength(pts, circuits * n_per_circuit + 1)
    return pts / np.linalg.norm(pts, axis=1).max()


def hopf_torus(eta: float = 0.70, circuits: int = 24, revolutions: float = 1.0, n_per_circuit: int = 150) -> np.ndarray:
    """The limit of the continuation: η fixed, φ sweeping `revolutions` × 360° over
    `circuits` fibres — a (circuits, revolutions) torus-knot conductor lying on the
    Hopf torus of that η.  This is the 'donut' the drift winding is heading toward."""

    n_hi = max(4097, circuits * 400 + 1)
    u = np.linspace(0.0, 1.0, n_hi)
    eta_arr = np.full_like(u, eta)
    phi = 2.0 * np.pi * revolutions * u
    tau = np.pi / 2.0 + 2.0 * np.pi * circuits * u
    pts = _hopf_points(eta_arr, phi, tau)
    pts = resample_arclength(pts, circuits * n_per_circuit + 1)
    return pts / np.linalg.norm(pts, axis=1).max()


# ---------------------------------------------------------------------------
# Multi-core assemblies
# ---------------------------------------------------------------------------


def rotation_matrix(axis, angle_rad: float) -> np.ndarray:
    a = unit(np.asarray(axis, dtype=float))
    c, s = math.cos(angle_rad), math.sin(angle_rad)
    x, y, z = a
    return np.array([[c + x * x * (1 - c), x * y * (1 - c) - z * s, x * z * (1 - c) + y * s],
                     [y * x * (1 - c) + z * s, c + y * y * (1 - c), y * z * (1 - c) - x * s],
                     [z * x * (1 - c) - y * s, z * y * (1 - c) + x * s, c + z * z * (1 - c)]])


ASSEMBLY_LAYOUTS = {
    # name: list of (offset direction, distance, rotation axis, rotation angle, current sign)
    "single": [((0, 0, 0), 0.0, (0, 0, 1), 0.0, +1)],
    "pair-axial-same": [((0, 0, 1), 0.5, (0, 0, 1), 0.0, +1), ((0, 0, -1), 0.5, (0, 0, 1), 0.0, +1)],
    "pair-axial-opposed": [((0, 0, 1), 0.5, (0, 0, 1), 0.0, +1), ((0, 0, -1), 0.5, (0, 0, 1), 0.0, -1)],
    "pair-facing-flipped": [((0, 0, 1), 0.5, (1, 0, 0), math.pi, +1), ((0, 0, -1), 0.5, (0, 0, 1), 0.0, +1)],
    "pair-side-by-side": [((1, 0, 0), 0.5, (0, 0, 1), 0.0, +1), ((-1, 0, 0), 0.5, (0, 0, 1), 0.0, +1)],
    "pair-side-crossed": [((1, 0, 0), 0.5, (0, 0, 1), 0.0, +1), ((-1, 0, 0), 0.5, (1, 0, 0), math.pi / 2, +1)],
    "triad-120": [((math.cos(k * 2 * math.pi / 3), math.sin(k * 2 * math.pi / 3), 0), 0.55, (0, 0, 1), k * 2 * math.pi / 3, +1) for k in range(3)],
    "triad-120-alternating": [((math.cos(k * 2 * math.pi / 3), math.sin(k * 2 * math.pi / 3), 0), 0.55, (0, 0, 1), k * 2 * math.pi / 3, (+1 if k % 2 == 0 else -1)) for k in range(3)],
    "tetra-4": [(d, 0.5, (0, 0, 1), 0.0, +1) for d in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))],
    "octa-6-cusp": [((1, 0, 0), 0.55, (0, 1, 0), math.pi / 2, +1), ((-1, 0, 0), 0.55, (0, 1, 0), -math.pi / 2, +1),
                    ((0, 1, 0), 0.55, (1, 0, 0), -math.pi / 2, -1), ((0, -1, 0), 0.55, (1, 0, 0), math.pi / 2, -1),
                    ((0, 0, 1), 0.55, (0, 0, 1), 0.0, +1), ((0, 0, -1), 0.55, (1, 0, 0), math.pi, +1)],
    "octa-6-aligned": [((1, 0, 0), 0.55, (0, 1, 0), math.pi / 2, +1), ((-1, 0, 0), 0.55, (0, 1, 0), -math.pi / 2, +1),
                       ((0, 1, 0), 0.55, (1, 0, 0), -math.pi / 2, +1), ((0, -1, 0), 0.55, (1, 0, 0), math.pi / 2, +1),
                       ((0, 0, 1), 0.55, (0, 0, 1), 0.0, +1), ((0, 0, -1), 0.55, (1, 0, 0), math.pi, +1)],
    "cube-8": [(d, 0.55, (0, 0, 1), 0.0, +1) for d in ((1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1), (-1, 1, 1), (-1, 1, -1), (-1, -1, 1), (-1, -1, -1))],
    "cube-8-checker": [(d, 0.55, (0, 0, 1), 0.0, (+1 if (d[0] * d[1] * d[2]) > 0 else -1)) for d in ((1, 1, 1), (1, 1, -1), (1, -1, 1), (1, -1, -1), (-1, 1, 1), (-1, 1, -1), (-1, -1, 1), (-1, -1, -1))],
}


def assembly(core_active: np.ndarray, layout: str, core_scale: float = 0.45) -> list:
    """Place scaled copies of a unit-ball core according to a layout.

    Returns a list of (points, current_sign) for each core, all inside the unit
    ball (core radius `core_scale`, centres at the layout distance).  Each core is
    a separate closed circuit; use `assembly_windings` to attach returns.
    """

    spec = ASSEMBLY_LAYOUTS[layout]
    cores = []
    for direction, dist, axis, angle, sign in spec:
        d = unit(np.asarray(direction, dtype=float)) if np.linalg.norm(direction) > 0 else np.zeros(3)
        Rm = rotation_matrix(axis, angle)
        pts = (core_active * (core_scale if len(spec) > 1 else 1.0)) @ Rm.T + dist * d
        cores.append((pts, sign))
    return cores


def assembly_windings(core_active: np.ndarray, layout: str, scale_m: float, wire_radius_m: float, current_A: float,
                      name: str, return_radius_factor: float = 3.8, core_scale: float = 0.45) -> list:
    """Closed Windings (one per core) for an assembly; the field is the sum over cores."""

    out = []
    for k, (pts, sign) in enumerate(assembly(core_active, layout, core_scale)):
        closed, n_active = close_with_return(pts * scale_m, return_radius_factor * scale_m)
        out.append(Winding(f"{name} [{layout} #{k}]", closed, sign * current_A, wire_radius_m, closed=True, active_count=n_active,
                           meta={"family": "assembly", "layout": layout, "core": name, "index": k, "sign": sign}))
    return out


class MultiWinding:
    """Several closed Windings treated as one source (fields add).  Quacks like a Winding for the field code."""

    def __init__(self, name: str, windings: list, meta: Optional[dict] = None):
        self.name = name
        self.windings = list(windings)
        self.meta = meta or {}
        self.wire_radius_m = min(w.wire_radius_m for w in windings)
        self.current_A = windings[0].current_A
        self.points = np.vstack([w.points for w in windings])
        self.active_count = None

    @property
    def active_points(self) -> np.ndarray:
        return np.vstack([w.active_points for w in self.windings])

    @property
    def length_m(self) -> float:
        return sum(w.length_m for w in self.windings)

    @property
    def active_length_m(self) -> float:
        return sum(w.active_length_m for w in self.windings)

    def clearance_m(self, skip: int = 12) -> float:
        own = min(w.clearance_m(skip) for w in self.windings)
        cross = float("inf")
        for i in range(len(self.windings)):
            tree = cKDTree(self.windings[i].active_points)
            for j in range(i + 1, len(self.windings)):
                cross = min(cross, float(tree.query(self.windings[j].active_points)[0].min()))   # every point (a touch at one point counts)
        return min(own, cross)

    def scaled(self, factor: float, name: Optional[str] = None):
        return MultiWinding(name or self.name, [w.scaled(factor) for w in self.windings], dict(self.meta))

    def with_current(self, current_A: float):
        """Set the physical current of the reference winding (the first) to |current_A| and scale every
        other winding by its relative amplitude, so relative currents are preserved (before v2.3 every
        winding was set to ±current_A, which silently equalised amplitudes in the Joule-power report)."""
        ref = abs(self.windings[0].current_A) if self.windings and self.windings[0].current_A != 0 else 1.0
        return MultiWinding(self.name, [w.with_current(current_A * w.current_A / ref) for w in self.windings], dict(self.meta))

    def with_wire_radius(self, a: float):
        return MultiWinding(self.name, [w.with_wire_radius(a) for w in self.windings], dict(self.meta))

    def resistance_ohm(self, resistivity_ohm_m: float = 1.68e-8) -> float:
        return sum(w.resistance_ohm(resistivity_ohm_m) for w in self.windings)

    def joule_power_W(self, resistivity_ohm_m: float = 1.68e-8) -> float:
        return sum(w.joule_power_W(resistivity_ohm_m) for w in self.windings)


# ---------------------------------------------------------------------------
# Hopf torus + mirror image: the "meshed torus"
# ---------------------------------------------------------------------------


def _weave_profile(u_cross: np.ndarray, sign: np.ndarray, u: np.ndarray, periodic: bool = False) -> np.ndarray:
    """Smooth ±1 profile along a curve that takes the value `sign` at each crossing
    parameter `u_cross` (cosine ramps between crossings).  periodic=True treats u as a
    closed parameter of period 1 (closed strands)."""

    order = np.argsort(u_cross)
    uc, s = u_cross[order], sign[order]
    if periodic:
        uu = np.r_[uc - 1.0, uc, uc + 1.0]
        ss = np.r_[s, s, s]
    else:
        uu = np.r_[uc[0] - (uc[1] - uc[0]), uc, uc[-1] + (uc[-1] - uc[-2])]
        ss = np.r_[-s[0], s, -s[-1]]
    i = np.clip(np.searchsorted(uu, u) - 1, 0, len(uu) - 2)
    t = np.clip((u - uu[i]) / (uu[i + 1] - uu[i]), 0.0, 1.0)
    return ss[i] + (ss[i + 1] - ss[i]) * (1.0 - np.cos(np.pi * t)) / 2.0


def hopf_mirror_pair(eta: float = 0.70, circuits: int = 24, revolutions: float = 1.0, mode: str = "woven",
                     delta_eta: float = 0.025, n_per_circuit: int = 150, tau0: float = 0.0, deform=None,
                     mirror_rotation: float = math.pi, info: Optional[dict] = None) -> list:
    """A Hopf-torus conductor and its mirror image (z → −z), fitted together on the
    same torus.  The mirror of a Hopf fibre is an anti-Hopf fibre, so the two
    families cross each other 2·(circuits² − revolutions²/4) times and make a
    rhombic mesh on the torus.

    mode="woven"  : each conductor is pushed alternately outward / inward (η ± δη)
                    at successive crossings, so the two interlock like a plain weave —
                    the "meshed torus".  Crossings are enumerated analytically:
                    with p = N + r/2, q = N − r/2 they sit at u = (j/p + k/q)/2 on the
                    original and v = (k/q − j/p)/2 on the mirror; the over/under
                    pattern (−1)^k alternates along both strands.
    mode="nested" : the original lies on the torus η + δη and the mirror on η − δη
                    (no interlocking; the two just nest).

    tau0 is the fibre phase at which each strand starts and ends.  Both ends of a
    strand sit at the same poloidal position, and that one missing turn leaves a
    φ-averaged radial field ∝ sin(tau0)/N at the midplane (1.2 % of B_tor for
    tau0 = π/2, N = 24); tau0 = 0 puts the ends on the midplane and cancels it, so
    that is the default (the session-5 sweep used π/2).

    deform: optional callable pts -> pts applied to both strands before the final
    normalisation (see `torus_deformation` — rotating-ellipse / helical-axis twist).

    mirror_rotation α: the mirror strand is B = R_z(α)·M·A.  With α = π (default) the
    weave sign (−1)^j is consistent with B being the *exact* rotated mirror image of the
    displaced strand A, so the pair's toroidal currents cancel identically and no
    stray poloidal field is left (the session-5 sweep used α = 0 with sign (−1)^k, where
    A and B carry different over/under profiles and a ~1 % residual B_z remains).

    Returns [points_A, points_B] in the unit ball.  With equal currents the
    poloidal components cancel and the pair acts like a toroidal current sheet;
    with opposed currents the toroidal components cancel and it acts like a
    toroidal-field coil (closed field lines inside the tube).

    info: optional dict, filled with the build details (p, q, normalisation scale, crossing parameters
    u along A / v along B, the over/under sign at each crossing, and whether the single-profile
    (exact-mirror) scheme was used).
    """

    N, r = circuits, revolutions
    if info is None:
        info = {}
    n_hi = max(4097, circuits * 400 + 1)
    u = np.linspace(0.0, 1.0, n_hi)
    alpha = mirror_rotation
    ca, sa = math.cos(alpha), math.sin(alpha)

    def mirror_of(P):   # R_z(α) · diag(1, 1, −1)
        return np.column_stack((ca * P[:, 0] - sa * P[:, 1], sa * P[:, 0] + ca * P[:, 1], -P[:, 2]))

    def curve(eta_arr, uu):
        return _hopf_points(eta_arr, 2.0 * np.pi * r * uu, tau0 + 2.0 * np.pi * N * uu)

    if mode == "nested":
        A = curve(np.full_like(u, eta + delta_eta), u)
        B = mirror_of(curve(np.full_like(u, eta - delta_eta), u))
    elif mode == "woven":
        p, q = N + r / 2.0, N - r / 2.0
        J = np.arange(-int(p) - 3, int(p) + 3)
        K = np.arange(-1, int(2 * q) + 3)
        jj, kk = np.meshgrid(J, K, indexing="ij")
        jj, kk = jj.ravel(), kk.ravel()
        # crossings of A(u) with R_α M A(v): a1(u) = a1(v) + α + 2πj,  a2(u) + a2(v) = π + 2πk
        uc = ((jj + alpha / (2 * np.pi)) / p + (kk + 0.5 - tau0 / np.pi) / q) / 2.0
        vc = ((kk + 0.5 - tau0 / np.pi) / q - (jj + alpha / (2 * np.pi)) / p) / 2.0
        closed_strand = abs(revolutions - round(revolutions)) < 1e-9 and int(round(revolutions)) % 2 == 0
        hi = 1.0 - 1e-12 if closed_strand else 1.0
        ok = (uc >= 0) & (uc <= hi) & (vc >= 0) & (vc <= hi)
        uc, vc, jj, kk = uc[ok], vc[ok], jj[ok], kk[ok]
        # B = R_α M A exactly: the crossing at u on A meets the mirror image of the crossing at v on A, so a
        # single over/under profile s(u) must satisfy s(u) = −s(v).  Alternating by rank along A does that
        # for (N even, revolutions = 2, α = π): every crossing pair has odd rank sum.  Check, else fall back.
        rank = np.empty(len(uc), int)
        order = np.argsort(uc)
        rank[order] = np.arange(len(uc))
        us = uc[order]
        partner = order[np.clip(np.searchsorted(us, vc), 0, len(us) - 1)]
        consistent = np.all(np.abs(uc[partner] - vc) < 1e-9) and np.all((rank + rank[partner]) % 2 == 1)
        if consistent:
            sign = (-1.0) ** rank
            A = curve(eta + delta_eta * _weave_profile(uc, sign, u, closed_strand), u)
            B = mirror_of(A)
        else:
            sign = (-1.0) ** kk
            A = curve(eta + delta_eta * _weave_profile(uc, sign, u, closed_strand), u)
            B = mirror_of(curve(eta + delta_eta * _weave_profile(vc, -sign, u, closed_strand), u))
        info.update({"crossing_u": uc.tolist(), "crossing_v": vc.tolist(), "crossing_j": jj.tolist(), "crossing_k": kk.tolist(),
                     "crossing_sign": [int(s) for s in sign], "exact_mirror": bool(consistent)})
    else:
        raise ValueError("mode must be 'woven' or 'nested'")
    if deform is not None:
        A, B = deform(A), deform(B)
    scale = 1.0 / max(np.linalg.norm(A, axis=1).max(), np.linalg.norm(B, axis=1).max())
    info.update({"p": float(N + r / 2.0), "q": float(N - r / 2.0), "strand_scale": float(scale), "tau0": float(tau0),
                 "mirror_rotation": float(alpha)})
    out = []
    for pts in (A, B):
        pts = resample_arclength(pts * scale, circuits * n_per_circuit + 1)
        out.append(pts)
    return out


def hopf_mirror_windings(scale_m: float, wire_radius_m: float, current_A: float, sense: int = 1,
                         name: str = "hopf mirror pair", return_radius_factor: float = 3.8, **kw) -> "MultiWinding":
    """Closed two-circuit MultiWinding for the mirror pair; sense=+1 same current,
    sense=-1 opposed (mirror carries −I)."""

    pair = hopf_mirror_pair(**kw)
    ws = []
    for k, (pts, sgn) in enumerate(zip(pair, (1, sense))):
        closed, n_active = close_with_return(pts * scale_m, return_radius_factor * scale_m)
        ws.append(Winding(f"{name} #{k}", closed, sgn * current_A, wire_radius_m, closed=True, active_count=n_active,
                          meta={"family": "hopf-mirror", "index": k, "sign": sgn, **kw}))
    return MultiWinding(name, ws, {"family": "hopf-mirror", "sense": sense, **kw})


# ---------------------------------------------------------------------------
# New geometries (session 5): baseball seam, yin-yang, picket fence, sphere winding, Hopf link
# ---------------------------------------------------------------------------


def baseball_seam(radius: float = 0.72, amplitude_deg: float = 40.0, turns: int = 6, radial_pitch: float = 0.035,
                  n_per_turn: int = 720) -> np.ndarray:
    """Baseball-seam ('yin-yang' half) coil: the closed curve on a sphere with
    latitude λ(t) = A·sin(2t), longitude t.  Wound `turns` times with the sphere
    radius stepping inward by `radial_pitch` each turn (nested seams), so the
    conductor never touches itself.  Classic minimum-B mirror coil.  Unit ball."""

    A = np.deg2rad(amplitude_deg)
    t = np.linspace(0.0, 2.0 * np.pi * turns, turns * n_per_turn + 1)
    lam = A * np.sin(2.0 * t)
    r = radius - radial_pitch * (t / (2.0 * np.pi))
    pts = np.column_stack((r * np.cos(lam) * np.cos(t), r * np.cos(lam) * np.sin(t), r * np.sin(lam)))
    return pts


def yin_yang(radius_outer: float = 0.76, radius_inner: float = 0.60, amplitude_deg: float = 40.0, turns: int = 4,
             radial_pitch: float = 0.03) -> list:
    """Two baseball coils at the same centre, the inner one rotated 90° about z so the
    lobes interleave (Livermore yin-yang).  Returns [outer_pts, inner_pts]; the
    two are separate circuits."""

    outer = baseball_seam(radius_outer, amplitude_deg, turns, radial_pitch)
    inner = baseball_seam(radius_inner, amplitude_deg, turns, radial_pitch) @ rotation_matrix([0, 0, 1], np.pi / 2).T
    return [outer, inner]


def picket_fence(rings: int = 5, radius: float = 0.72, n: int = 720) -> list:
    """Stack of coaxial rings on the sphere of `radius`, equally spaced in z, with
    alternating current direction (picket-fence / ring-cusp).  Returns
    [(points, sign), ...] — one circuit per ring."""

    z = np.linspace(-0.62, 0.62, rings) * radius
    out = []
    for k, zk in enumerate(z):
        rk = math.sqrt(max(radius**2 - zk**2, 1e-6))
        pts = circular_loop(rk, n, center=(0.0, 0.0, float(zk)))
        out.append((pts, (-1) ** k))
    return out


def sphere_winding(radius: float = 0.76, turns: int = 24, n_per_turn: int = 90) -> np.ndarray:
    """Spherical solenoid: one conductor spiralling from pole to pole on the sphere
    with turns equally spaced in z.  The surface current then goes as sin θ, which
    gives an exactly uniform interior field B = μ0·N·I/(3R) and a pure dipole
    outside — a clean analytic check as well as a 'cool' shape.  Unit ball."""

    t = np.linspace(0.0, 1.0, turns * n_per_turn + 1)
    z = radius * (0.94 * (2.0 * t - 1.0))
    rho = np.sqrt(np.maximum(radius**2 - z**2, 0.0))
    ang = 2.0 * np.pi * turns * t
    return np.column_stack((rho * np.cos(ang), rho * np.sin(ang), z))


def cone_helix(radius_base: float = 0.32, radius_tip: float = 0.06, height: float = 0.5, turns: int = 8,
               n_per_turn: int = 90) -> np.ndarray:
    """'Tornado' coil: a conical helix along z, base (radius_base) at z = -height/2, tip
    (radius_tip) at z = +height/2, `turns` turns equally spaced in z.  The on-axis field of
    each turn is ∝ 1/R at its own plane, so |B| rises toward the tip: a converging
    (mirror-like) field, the closest static magnetostatics gets to a 'funnel'.  Unit ball;
    orient/place it with the design transform (reflect z to point the tip down)."""

    t = np.linspace(0.0, 1.0, int(turns) * n_per_turn + 1)
    r = radius_base + (radius_tip - radius_base) * t
    z = -0.5 * height + height * t
    ang = 2.0 * np.pi * turns * t
    return np.column_stack((r * np.cos(ang), r * np.sin(ang), z))


def torus_helical_windings(R0: float = 0.60, r: float = 0.40, l: int = 2, periods: int = 4, phase_deg: float = 0.0,
                           n_per_turn: int = 720) -> list:
    """Classical stellarator helical windings: 2l closed conductors on the torus (R0, r) following
    θ = (n/l)·φ + kπ/l (k = 0 … 2l−1), adjacent conductors carrying opposite currents, n field periods.
    Returns [(points, sign), ...]; each conductor closes on itself after l/gcd(n, l) toroidal turns, so
    there are 2·gcd(n, l) distinct closed conductors (4 for l = 2, n even; 2 for n odd).
    With a toroidal field inside (e.g. the meshed torus) this is the l = 2 stellarator: rotational
    transform from the external windings alone, no weave modulation needed.  Unit ball if R0 + r ≤ 1."""

    from math import gcd
    g = gcd(abs(int(periods)), l)
    m = l // g                                     # toroidal turns until a conductor closes
    out = []
    for k in range(2 * g):                         # distinct closed conductors (2l / m); k and k + 2g coincide otherwise
        phi = np.linspace(0.0, 2.0 * np.pi * m, m * n_per_turn + 1)
        theta = (periods / l) * phi + k * np.pi / l + np.deg2rad(phase_deg)
        pts = np.column_stack(((R0 + r * np.cos(theta)) * np.cos(phi), (R0 + r * np.cos(theta)) * np.sin(phi), r * np.sin(theta)))
        pts[-1] = pts[0]
        out.append((pts, (-1) ** k))
    return out


def hopf_link(radius: float = 0.62, n: int = 720) -> list:
    """Two linked circular loops (Hopf link): loop A in the xy-plane centred at
    (−R/2, 0, 0); loop B in the xz-plane centred at (+R/2, 0, 0)."""

    a = circular_loop(radius, n, center=(-radius / 2, 0.0, 0.0), axis="z")
    b = circular_loop(radius, n, center=(radius / 2, 0.0, 0.0), axis="y")
    return [a, b]


def multi_windings(parts: list, scale_m: float, wire_radius_m: float, current_A: float, name: str, family: str,
                   return_radius_factor: float = 3.8, meta: Optional[dict] = None) -> "MultiWinding":
    """Close each (points, sign) part with a remote return and bundle them."""

    ws = []
    for k, part in enumerate(parts):
        pts, sgn = part if isinstance(part, tuple) else (part, 1)
        closed, n_active = close_with_return(np.asarray(pts) * scale_m, return_radius_factor * scale_m)
        ws.append(Winding(f"{name} #{k}", closed, sgn * current_A, wire_radius_m, closed=True, active_count=n_active,
                          meta={"family": family, "index": k, "sign": sgn}))
    return MultiWinding(name, ws, {"family": family, **(meta or {})})


def torus_deformation(R0: float, r_tube: float, ellipticity: float = 0.0, periods: int = 5, helical_axis: float = 0.0,
                      triangularity: float = 0.0):
    """Map that twists a round torus (major radius R0, tube radius r_tube, axis z) into a
    stellarator-like winding surface:

      r(θ, ζ) = r · [1 + ε cos(2θ − nζ) + δ cos(3θ − nζ)]   rotating ellipse (l = 2) and
                                                             triangle (l = 3), n periods
      axis    = (R0 + h·r_tube cos nζ, h·r_tube sin nζ)      helical magnetic-axis excursion

    Returns a callable pts -> pts.  With ε = h = δ = 0 it is the identity.
    """

    n = periods

    def f(pts: np.ndarray) -> np.ndarray:
        x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
        rho = np.hypot(x, y)
        zeta = np.arctan2(y, x)
        dr = rho - R0
        r = np.hypot(dr, z)
        th = np.arctan2(z, dr)
        r2 = r * (1.0 + ellipticity * np.cos(2 * th - n * zeta) + triangularity * np.cos(3 * th - n * zeta))
        Rax = R0 + helical_axis * r_tube * np.cos(n * zeta)
        zax = helical_axis * r_tube * np.sin(n * zeta)
        rho2 = Rax + r2 * np.cos(th)
        z2 = zax + r2 * np.sin(th)
        return np.column_stack((rho2 * np.cos(zeta), rho2 * np.sin(zeta), z2))

    return f


def hopf_torus_frame(eta: float = 0.70) -> tuple:
    """(R0, r_tube) of the round torus that hopf_torus / hopf_mirror_pair(eta) lie on, in the
    normalisation used by those generators (outer edge at radius 1)."""

    u = np.linspace(0.0, 1.0, 4001)
    pts = _hopf_points(np.full_like(u, eta), np.zeros_like(u), 2 * np.pi * u)
    rho = np.hypot(pts[:, 0], pts[:, 1])
    s = 1.0 / np.linalg.norm(pts, axis=1).max()
    return float(0.5 * (rho.min() + rho.max()) * s), float(0.5 * (rho.max() - rho.min()) * s)


def twisted_mirror_pair(eta: float = 0.70, circuits: int = 24, ellipticity: float = 0.2, periods: int = 5,
                        helical_axis: float = 0.0, triangularity: float = 0.0, **kw) -> list:
    """hopf_mirror_pair on a twisted winding surface (see torus_deformation)."""

    R0, rt = hopf_torus_frame(eta)
    # the generator normalises after deformation; deform in the pre-normalised frame, so rescale R0, rt
    u = np.linspace(0.0, 1.0, 4001)
    raw = _hopf_points(np.full_like(u, eta), np.zeros_like(u), 2 * np.pi * u)
    s = np.linalg.norm(raw, axis=1).max()
    f = torus_deformation(R0 * s, rt * s, ellipticity, periods, helical_axis, triangularity)
    return hopf_mirror_pair(eta, circuits, mode=kw.pop("mode", "woven"), deform=f, **kw)


def _modulated_strand(eta: float, p: float, q: float, eps: float, n: int, phase: float, handed: int, n_hi: int):
    """Strand on the Hopf torus η whose turns are the level sets of
        F(a1, a2) = a2 + (eps/2)·sin(2·a2 − handed·n·a1 + phase) − (q/p)·a1,
    so the turn density (∝ surface current) is ∂F/∂a2 = 1 + eps·cos(2 a2 − handed n a1 + phase):
    an l = 2, n-period helical modulation.  (A modulation of the local *rate* da2/da1
    averages out over each toroidal transit and leaves the turn spacing unchanged — the
    displacement form is what actually moves current around the tube.)  Returns (a1, a2)
    on a1 ∈ [0, 2πp]; the strand closes on itself after 2πp since sin(...) is periodic."""

    a1 = np.linspace(0.0, 2 * np.pi * p, n_hi)
    a2 = (q / p) * a1
    # fixed-point iteration for a2 = (q/p) a1 − (eps/2) sin(2 a2 − handed n a1 + phase): a contraction with
    # factor |eps| < 1, iterated to convergence (12 sweeps left a 0.7^12 ≈ 1 % residual at eps = 0.7)
    for _ in range(2000):
        a2_new = (q / p) * a1 - 0.5 * eps * np.sin(2 * a2 - handed * n * a1 + phase)
        done = np.max(np.abs(a2_new - a2)) < 1e-13
        a2 = a2_new
        if done:
            break
    return a1, a2


def _eta_rate(a1: np.ndarray, a2: np.ndarray, eta: float, pts, h: float = 1e-4) -> np.ndarray:
    """|∂p/∂η| along a strand: how far a unit change of η moves the conductor off the winding surface."""
    return np.linalg.norm(pts(a1, a2, np.full_like(a1, eta + h)) - pts(a1, a2, np.full_like(a1, eta - h)), axis=1) / (2 * h)


def _runs(mask: np.ndarray) -> list:
    """Maximal runs of True in a periodic boolean array → list of (start, stop) index pairs, stop exclusive;
    a run that wraps is returned with stop > len(mask) (indices taken mod len)."""
    n = len(mask)
    if mask.all():
        return [(0, n)]
    if not mask.any():
        return []
    # rotate so that index 0 is free, then runs never wrap
    k = int(np.argmin(mask))             # first False
    m = np.roll(mask, -k)
    d = np.diff(np.r_[False, m, False].astype(int))
    starts, stops = np.where(d == 1)[0], np.where(d == -1)[0]
    return [(int(a + k), int(b + k)) for a, b in zip(starts, stops)]


def _zone_weave(A0: np.ndarray, B0: np.ndarray, rateA: np.ndarray, rateB: np.ndarray, delta_eta: float, u: np.ndarray,
                uc: np.ndarray, vc: np.ndarray, sign: np.ndarray, margin: float = 1.15) -> tuple:
    """Over/under profiles for two strands on the same surface that keep them apart everywhere.

    A crossing-based weave (±1 at each crossing, cosine ramps between) fails wherever the two strands
    run close together *in the surface* without crossing, or cross several times within a short
    stretch: a ramp through zero there puts both conductors at η and they touch.  So the profile is
    built from *contact zones* instead: wherever the in-surface distance g between the strands is
    below the normal separation h = 2 δη |∂p/∂η| the two profiles are frozen at opposite full values,
    and ramps are allowed only in the free stretches between zones (spanning the whole stretch).
    Zones of A and B that face each other are joined into one component (union-find), so every
    contact has A over and B under (or vice versa); the component's sign is the plain-weave
    (checkerboard) sign of its earliest crossing, which keeps the alternation wherever the strands
    are far enough apart to allow it.  Returns (sA, sB, info)."""

    n = len(u)
    tB, tA = cKDTree(B0), cKDTree(A0)
    gA = tB.query(A0)[0]
    gB = tA.query(B0)[0]
    hA = 2 * abs(delta_eta) * rateA
    hB = 2 * abs(delta_eta) * rateB
    zonesA = _runs(gA < margin * hA)
    zonesB = _runs(gB < margin * hB)
    zidA = -np.ones(n, int)
    for k, (a, b) in enumerate(zonesA):
        zidA[np.arange(a, b) % n] = k
    zidB = -np.ones(n, int)
    for k, (a, b) in enumerate(zonesB):
        zidB[np.arange(a, b) % n] = k
    nA, nB = len(zonesA), len(zonesB)
    parent = list(range(nA + nB))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[max(rx, ry)] = min(rx, ry)

    # every B point within reach of an A zone point faces it: join their zones (and the symmetric pass)
    inA = np.where(zidA >= 0)[0]
    for i, nb in zip(inA, tB.query_ball_point(A0[inA], hA[inA])):
        for j in nb:
            if zidB[j] >= 0:
                union(zidA[i], nA + zidB[j])
    inB = np.where(zidB >= 0)[0]
    for j, nb in zip(inB, tA.query_ball_point(B0[inB], hB[inB])):
        for i in nb:
            if zidA[i] >= 0:
                union(zidA[i], nA + zidB[j])
    # component signs: the checkerboard sign of the earliest crossing in the component
    comp_sign: dict = {}
    order = np.argsort(uc)
    for k in order:
        ia = min(n - 1, int(round(uc[k] * (n - 1))))
        za = zidA[ia]
        if za < 0:                       # a crossing always sits in a zone (g = 0 there); guard for rounding
            continue
        c = find(za)
        comp_sign.setdefault(c, float(sign[k]))
    crossing_sign = np.array([comp_sign.get(find(zidA[min(n - 1, int(round(uu * (n - 1))))]), s) if zidA[min(n - 1, int(round(uu * (n - 1))))] >= 0 else s
                              for uu, s in zip(uc, sign)])
    defects = int(np.sum(crossing_sign != sign))

    def profile(zones, zid, comp_of_zone, own_sign_factor):
        # zone signs (components without a crossing inherit the previous zone's sign along the strand)
        zs = np.zeros(len(zones))
        prev = 1.0
        for k in range(len(zones)):
            c = comp_of_zone(k)
            if c in comp_sign:
                zs[k] = comp_sign[c]
            else:
                zs[k] = prev
                comp_sign[c] = prev
            prev = zs[k]
        out = np.ones(n)
        if not zones:
            return out * own_sign_factor
        if len(zones) == 1 and zones[0][1] - zones[0][0] >= n:
            return out * zs[0] * own_sign_factor
        for k, (a, b) in enumerate(zones):
            out[np.arange(a, b) % n] = zs[k]
            # free stretch from this zone's end to the next zone's start
            a2, _ = zones[(k + 1) % len(zones)]
            if a2 <= b:
                a2 += n
            s0, s1 = zs[k], zs[(k + 1) % len(zones)]
            idx = np.arange(b, a2)
            if len(idx):
                t = (idx - b + 1.0) / (len(idx) + 1.0)
                out[idx % n] = s0 + (s1 - s0) * (1.0 - np.cos(np.pi * t)) / 2.0 if s0 != s1 else s0
        return out * own_sign_factor

    sA = profile(zonesA, zidA, lambda k: find(k), 1.0)
    sB = profile(zonesB, zidB, lambda k: find(nA + k), -1.0)
    info = {"defects": defects, "zone_fraction": float(np.mean(zidA >= 0)), "components": int(len({find(x) for x in range(nA + nB)})),
            "crossing_sign": crossing_sign, "min_distance": float("nan")}
    return sA, sB, info


def hopf_helical_pair(eta: float = 0.70, circuits: int = 24, eps: float = 0.2, periods: int = 5, delta_eta: float = 0.025,
                      mirror_rotation: float = math.pi, n_per_circuit: int = 150, n_hi: int = 20000, deform=None,
                      helical_mode: str = "toroidal", return_crossings: bool = False) -> tuple:
    """The meshed torus made chiral: both strand families carry the same-handed helical
    density modulation 1 + eps·cos(2θ − nζ) (an l = 2, n-period stellarator surface
    current).  Family B is the rotated mirror image of a strand whose modulation was
    pre-mirrored, so after mirroring both families have the same handedness; their
    toroidal currents still cancel locally (equal densities, opposite currents) and
    their poloidal currents add with the helical modulation.  Crossings are found
    numerically and woven by rank along A.  Returns ([A, B], info).

    helical_mode="poloidal": both families modulated alike → the *poloidal* current (hence the
    toroidal field) carries the cos(2θ − nζ) pattern: a mirror-like ripple, little transform.
    helical_mode="toroidal": the families are modulated oppositely (A: 1+ε cos, B: 1−ε cos) →
    the poloidal current stays uniform and the *toroidal* current acquires the pattern
    2ε cos(2θ − nζ): the l = 2 helical-winding current of a classical stellarator, which
    is what produces rotational transform."""

    N = circuits
    p, q = N + 1.0, N - 1.0        # revolutions = 2: closed (N+1, N−1) torus knots
    alpha = mirror_rotation
    ca, sa = math.cos(alpha), math.sin(alpha)

    def mirror_of(P):
        return np.column_stack((ca * P[:, 0] - sa * P[:, 1], sa * P[:, 0] + ca * P[:, 1], -P[:, 2]))

    a1A, a2A = _modulated_strand(eta, p, q, eps, periods, 0.0, +1, n_hi)
    # B0 gets the opposite handedness (and the phase shift nα) so that R_α M B0 has the same modulation as A
    eps_B = eps if helical_mode == "poloidal" else -eps
    a1B, a2B = _modulated_strand(eta, p, q, eps_B, periods, -periods * alpha, -1, n_hi)

    def pts(a1, a2, eta_arr):
        return _hopf_points(eta_arr, a1 - a2, 0.5 * (a1 + a2))

    # --- crossings, exactly, in the (a1, a2) torus coordinates ---
    # B after mirroring sits at a1 = a1B + α, a2 = π − a2B.  A(s) meets B(t) when s ≡ t + α (mod 2π) and
    # a2A(s) + a2B(t) ≡ π (mod 2π).  Both strands are graphs over their own a1, so scan the 2p+1 branches.
    L = 2 * np.pi * p
    # closed strands: a2B(t + L) = a2B(t) + 2πq, so extend B periodically and let t run outside [0, L]
    a1B_ext = np.concatenate((a1B - L, a1B[1:], a1B[1:] + L))
    a2B_ext = np.concatenate((a2B - 2 * np.pi * q, a2B[1:], a2B[1:] + 2 * np.pi * q))
    # scan A one sample past each end (periodic extension) so a crossing sitting exactly on u = 0
    # is caught from whichever side floating point puts it — a missed one leaves an odd count and
    # a weave defect at the wrap-around
    a1A_s = np.concatenate(([a1A[-2] - L], a1A, [a1A[1] + L]))
    a2A_s = np.concatenate(([a2A[-2] - 2 * np.pi * q], a2A, [a2A[1] + 2 * np.pi * q]))
    uc_list, vc_list = [], []
    # branches m and m + p are the same crossing (a2B(t + 2πp) = a2B(t) + 2πq), so scan exactly the p
    # distinct branches; their t stays inside the extended B table (no clamped interpolation)
    for m in range(int(round(p))):
        t = a1A_s - alpha + 2 * np.pi * m
        a2Bt = np.interp(t, a1B_ext, a2B_ext)
        D = np.mod(a2A_s + a2Bt, 2 * np.pi) - np.pi
        ok = (np.sign(D[:-1]) != np.sign(D[1:])) & (np.abs(D[1:] - D[:-1]) < np.pi)
        i = np.where(ok)[0]
        f = D[i] / (D[i] - D[i + 1])
        s_cross = a1A_s[i] + f * (a1A_s[i + 1] - a1A_s[i])
        uc_list.append(s_cross / L)
        vc_list.append((s_cross - alpha + 2 * np.pi * m) / L)
    uc = np.concatenate(uc_list)
    vc = np.concatenate(vc_list)
    uc = np.mod(uc, 1.0)
    uc[uc > 1 - 1e-9] = 0.0            # a crossing sitting on u = 0 ranks first, never last
    vc = np.mod(vc, 1.0)
    # drop duplicates (the same crossing reached from two branches / both ends of the scan): merge
    # anything closer than 1e-6 in u (crossings are ~1/(2N²) apart), including across the wrap
    order = np.argsort(uc)
    uc, vc = uc[order], vc[order]
    keep = np.r_[True, np.diff(uc) > 1e-6]
    uc, vc = uc[keep], vc[keep]
    if len(uc) > 1 and (uc[0] + 1.0 - uc[-1]) < 1e-6:
        uc, vc = uc[:-1], vc[:-1]
    if len(uc) % 2:
        raise RuntimeError(f"odd crossing count {len(uc)}: weave parity broken (raise n_hi)")
    rank = np.empty(len(uc), int)
    rank[np.argsort(uc)] = np.arange(len(uc))
    sign = (-1.0) ** rank                      # plain-weave checkerboard: alternate by rank along A
    u = np.linspace(0.0, 1.0, n_hi)
    sA, sB, weave = _zone_weave(pts(a1A, a2A, np.full_like(a1A, eta)), mirror_of(pts(a1B, a2B, np.full_like(a1B, eta))),
                                _eta_rate(a1A, a2A, eta, pts), _eta_rate(a1B, a2B, eta, pts), delta_eta, u, uc, vc, sign)
    # alternation along B (diagnostic): plain weave if ~1
    ob = np.argsort(vc)
    sgn_used = weave["crossing_sign"]
    alt_B = float(np.mean(sgn_used[ob][1:] * sgn_used[ob][:-1] < 0)) if len(ob) > 1 else float("nan")
    A = pts(a1A, a2A, eta + delta_eta * sA)
    B = mirror_of(pts(a1B, a2B, eta + delta_eta * sB))
    weave["min_distance"] = float(cKDTree(A).query(B)[0].min())      # exact, on the fine curves
    if deform is not None:
        A, B = deform(A), deform(B)
    scale = 1.0 / max(np.linalg.norm(A, axis=1).max(), np.linalg.norm(B, axis=1).max())
    out = [resample_arclength(P * scale, N * n_per_circuit + 1) for P in (A, B)]
    for P in out:
        P[-1] = P[0]
    zf = weave["zone_fraction"]
    info = {"crossings": int(len(uc)), "alternation_B": alt_B, "weave_defects": weave["defects"], "zone_fraction": zf,
            "weave": "nested" if zf > 0.99 else ("woven" if zf < 0.5 else "mixed"),
            "components": weave["components"], "min_distance_fine": float(weave["min_distance"] * scale),
            "closure_check": float((a2A[-1] - a2A[0]) / (2 * np.pi * q))}
    if return_crossings:      # parameters of every crossing (u along A, v along B, over/under sign used) — plain lists
        info["crossing_u"] = uc.tolist()
        info["crossing_v"] = vc.tolist()
        info["crossing_sign"] = [int(s) for s in weave["crossing_sign"]]
        info["strand_scale"] = float(scale)
    return out, info
