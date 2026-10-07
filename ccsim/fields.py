"""Electromagnetic fields of closed-wire windings — the magnetostatic Maxwell
solution, its vector potential, Faraday induction, and the checks.

Fidelity statement
------------------
* **B** is the exact solution of ∇×B = μ0 J, ∇·B = 0 for a filamentary
  polyline (Hanson & Hirshman 2002 compact segment formula), with the
  interior of a round conductor of radius a treated as a uniform current
  (|B| ∝ ρ/a² inside).  No softening parameter is used outside the conductor.
* **A** is the exact vector potential of the same filaments (Coulomb gauge,
  regularised only inside the conductor), so B = ∇×A holds and the induced
  electric field of a current ramp, E_ind = −∂A/∂t = −(dI/dt)·A/I, satisfies
  Faraday's law ∇×E = −∂B/∂t identically.
* Displacement current is neglected (quasi-static: L/λ ≪ 1 for any ramp slower
  than nanoseconds at metre scale), and the plasma's own currents and space
  charge are *not* included — this is a vacuum-field + test-particle stack.
  Electrostatic bias fields can be superposed from point charges / a uniform E.
* Every claim above is checked numerically by `maxwell_checks` (∇·B, ∇×B
  outside the conductor, ∇×A − B, Ampère loop integral = μ0 I) and by
  `analytic_references` against closed-form loop / solenoid / wire fields.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Callable, Dict, Optional, Sequence, Tuple

import numpy as np
from scipy.special import ellipe, ellipk

from .constants import MU0
from .geometry import Winding, nearest_distance


# ---------------------------------------------------------------------------
# Exact segment fields
# ---------------------------------------------------------------------------


def biot_savart(points: np.ndarray, winding: Winding, chunk: int = 256, current_A: Optional[float] = None) -> np.ndarray:
    """Exact B (tesla) of a polyline current at `points` (M, 3) in metres."""

    if hasattr(winding, "windings"):   # MultiWinding: superpose
        parts = winding.windings
        return sum(biot_savart(points, w, chunk, None if current_A is None else current_A * (1 if w.current_A >= 0 else -1)) for w in parts)
    pts = np.atleast_2d(np.asarray(points, dtype=float))
    I = winding.current_A if current_A is None else current_A
    p1, dl = winding.segments()
    p2 = p1 + dl
    L = np.linalg.norm(dl, axis=1)
    keep = L > 1e-15
    p1, p2, dl, L = p1[keep], p2[keep], dl[keep], L[keep]
    that = dl / L[:, None]
    a_wire = winding.wire_radius_m
    out = np.zeros_like(pts)
    pref = MU0 * I / (4.0 * math.pi)
    tx, ty, tz = that[:, 0][None, :], that[:, 1][None, :], that[:, 2][None, :]
    for i in range(0, len(pts), chunk):
        r = pts[i : i + chunk]
        Rix = r[:, 0:1] - p1[None, :, 0]
        Riy = r[:, 1:2] - p1[None, :, 1]
        Riz = r[:, 2:3] - p1[None, :, 2]
        Rfx = r[:, 0:1] - p2[None, :, 0]
        Rfy = r[:, 1:2] - p2[None, :, 1]
        Rfz = r[:, 2:3] - p2[None, :, 2]
        ri = np.sqrt(Rix * Rix + Riy * Riy + Riz * Riz)
        rf = np.sqrt(Rfx * Rfx + Rfy * Rfy + Rfz * Rfz)
        dot = Rix * Rfx + Riy * Rfy + Riz * Rfz
        rirf = ri * rf
        denom = rirf * (rirf + dot)
        s0 = Rix * tx + Riy * ty + Riz * tz
        rho2 = np.maximum(ri * ri - s0 * s0, 0.0)
        inside = (rho2 < a_wire * a_wire) & (s0 >= -a_wire) & (s0 <= L[None, :] + a_wire)
        factor = np.where(inside, rho2 / (a_wire * a_wire), 1.0)
        safe = np.where(denom > 1e-300, denom, 1e-300)
        coef = np.where(denom > 1e-300, (ri + rf) / safe, 0.0) * factor
        cx = Riy * Rfz - Riz * Rfy
        cy = Riz * Rfx - Rix * Rfz
        cz = Rix * Rfy - Riy * Rfx
        out[i : i + chunk, 0] = pref * np.sum(cx * coef, axis=1)
        out[i : i + chunk, 1] = pref * np.sum(cy * coef, axis=1)
        out[i : i + chunk, 2] = pref * np.sum(cz * coef, axis=1)
    return out


def vector_potential(points: np.ndarray, winding: Winding, chunk: int = 256, current_A: Optional[float] = None) -> np.ndarray:
    """Exact Coulomb-gauge A (T·m) of the polyline current, regularised inside the conductor."""

    if hasattr(winding, "windings"):
        return sum(vector_potential(points, w, chunk, None if current_A is None else current_A * (1 if w.current_A >= 0 else -1)) for w in winding.windings)
    pts = np.atleast_2d(np.asarray(points, dtype=float))
    I = winding.current_A if current_A is None else current_A
    p1, dl = winding.segments()
    L = np.linalg.norm(dl, axis=1)
    keep = L > 1e-15
    p1, dl, L = p1[keep], dl[keep], L[keep]
    that = dl / L[:, None]
    a2 = winding.wire_radius_m**2
    out = np.zeros_like(pts)
    pref = MU0 * I / (4.0 * math.pi)
    for i in range(0, len(pts), chunk):
        r = pts[i : i + chunk]
        Ri = r[:, None, :] - p1[None, :, :]
        s0 = np.sum(Ri * that[None, :, :], axis=2)
        rho2 = np.maximum(np.sum(Ri * Ri, axis=2) - s0 * s0, 0.0)
        rho = np.sqrt(np.maximum(rho2, a2))   # exact outside the conductor, capped inside it
        term = np.arcsinh((L[None, :] - s0) / rho) - np.arcsinh(-s0 / rho)
        out[i : i + chunk] = pref * np.sum(term[:, :, None] * that[None, :, :], axis=1)
    return out


def induced_electric_field(points: np.ndarray, winding: Winding, dI_dt_A_per_s: float) -> np.ndarray:
    """E_ind = −∂A/∂t for a current ramp dI/dt (V/m).  Exactly Faraday-consistent."""

    return -dI_dt_A_per_s * vector_potential(points, winding, current_A=1.0)


# ---------------------------------------------------------------------------
# Electrostatic superposition (optional bias)
# ---------------------------------------------------------------------------


@dataclass
class Electrostatics:
    """Simple superposition: uniform field + point charges (Coulomb).  Enough for
    a biased grid / electrode approximation; not a Poisson solve."""

    uniform_E_V_per_m: Tuple[float, float, float] = (0.0, 0.0, 0.0)
    point_charges: Sequence[Tuple[Tuple[float, float, float], float]] = ()   # ((x,y,z), q [C])
    softening_m: float = 1e-3

    def field(self, points: np.ndarray) -> np.ndarray:
        pts = np.atleast_2d(np.asarray(points, dtype=float))
        E = np.tile(np.asarray(self.uniform_E_V_per_m, dtype=float), (len(pts), 1))
        k = 1.0 / (4.0 * math.pi * 8.8541878128e-12)
        for pos, q in self.point_charges:
            r = pts - np.asarray(pos, dtype=float)
            d2 = np.sum(r * r, axis=1) + self.softening_m**2
            E += k * q * r / d2[:, None] ** 1.5
        return E


# ---------------------------------------------------------------------------
# Grid cache with trilinear interpolation
# ---------------------------------------------------------------------------


class FieldGrid:
    """Uniform Cartesian grid caching B (and optionally A) of a winding.

    Trilinear interpolation is used by the particle pusher; the grid error is
    measured explicitly against direct evaluation (`interpolation_error`) so it
    can be reported rather than assumed.
    """

    def __init__(self, winding: Winding, half_extent_m: float, n: int = 41, with_A: bool = False, chunk: int = 256):
        if n < 5:
            raise ValueError("grid needs at least 5 points per axis")
        self.winding = winding
        self.half = float(half_extent_m)
        self.n = int(n)
        self.axis = np.linspace(-self.half, self.half, self.n)
        self.dx = float(self.axis[1] - self.axis[0])
        xx, yy, zz = np.meshgrid(self.axis, self.axis, self.axis, indexing="ij")
        self.points = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
        try:   # numba kernel when available (identical to biot_savart to ~1e-15, ~25× faster)
            from .fastfield import HAVE_NUMBA, SegmentSource
            fast = HAVE_NUMBA
        except Exception:   # pragma: no cover
            fast = False
        if fast:
            self.B = SegmentSource(winding).B(self.points).reshape(self.n, self.n, self.n, 3)
        else:
            self.B = biot_savart(self.points, winding, chunk=chunk).reshape(self.n, self.n, self.n, 3)
        self.A = vector_potential(self.points, winding, chunk=chunk).reshape(self.n, self.n, self.n, 3) if with_A else None
        self.Bmag = np.linalg.norm(self.B, axis=3)
        try:
            from scipy.spatial import cKDTree
            # every conductor point (v2.3): the earlier [::2] subsample made the distance grid — and through the
            # ROI mask the 2 T normalisation — depend on the concatenation order of the windings
            self.wire_distance = cKDTree(winding.points).query(self.points)[0].reshape(self.n, self.n, self.n)
        except Exception:   # pragma: no cover
            self.wire_distance = nearest_distance(self.points, winding.points).reshape(self.n, self.n, self.n)

    # -- interpolation --------------------------------------------------------
    def _weights(self, pts: np.ndarray):
        f = (pts + self.half) / self.dx
        i0 = np.clip(np.floor(f).astype(int), 0, self.n - 2)
        w = np.clip(f - i0, 0.0, 1.0)
        return i0, w

    def interp(self, grid: np.ndarray, pts: np.ndarray) -> np.ndarray:
        pts = np.atleast_2d(pts)
        i0, w = self._weights(pts)
        ix, iy, iz = i0[:, 0], i0[:, 1], i0[:, 2]
        wx, wy, wz = w[:, 0], w[:, 1], w[:, 2]
        vec = grid.ndim == 4
        out = np.zeros((len(pts), 3) if vec else len(pts))
        for dx in (0, 1):
            ax = wx if dx else 1.0 - wx
            for dy in (0, 1):
                ay = wy if dy else 1.0 - wy
                for dz in (0, 1):
                    az = wz if dz else 1.0 - wz
                    wt = ax * ay * az
                    vals = grid[ix + dx, iy + dy, iz + dz]
                    out += wt[:, None] * vals if vec else wt * vals
        return out

    def B_at(self, pts: np.ndarray) -> np.ndarray:
        return self.interp(self.B, pts)

    def A_at(self, pts: np.ndarray) -> np.ndarray:
        if self.A is None:
            raise ValueError("grid was built without A")
        return self.interp(self.A, pts)

    def Bmag_at(self, pts: np.ndarray) -> np.ndarray:
        return self.interp(self.Bmag, pts)

    def wire_distance_at(self, pts: np.ndarray) -> np.ndarray:
        return self.interp(self.wire_distance, pts)

    def inside(self, pts: np.ndarray) -> np.ndarray:
        return np.all(np.abs(pts) < self.half, axis=1)

    # -- checks ----------------------------------------------------------------
    def interpolation_error(self, n_samples: int = 200, seed: int = 0, min_wire_distance: Optional[float] = None) -> Dict[str, float]:
        rng = np.random.default_rng(seed)
        pts = rng.uniform(-0.8 * self.half, 0.8 * self.half, size=(n_samples * 3, 3))
        d = self.wire_distance_at(pts)
        thresh = (3.0 * self.dx) if min_wire_distance is None else min_wire_distance
        pts = pts[d > thresh][:n_samples]
        direct = biot_savart(pts, self.winding)
        interp = self.B_at(pts)
        rel = np.linalg.norm(interp - direct, axis=1) / np.maximum(np.linalg.norm(direct, axis=1), 1e-300)
        return {"n": int(len(pts)), "rel_p50": float(np.median(rel)), "rel_p95": float(np.percentile(rel, 95)), "rel_max": float(rel.max())}


# ---------------------------------------------------------------------------
# Maxwell checks
# ---------------------------------------------------------------------------


def _central_diff_field(func: Callable[[np.ndarray], np.ndarray], pts: np.ndarray, h: float):
    """Return Jacobian J[m, i, j] = dF_i/dx_j by central differences."""

    pts = np.atleast_2d(pts)
    J = np.zeros((len(pts), 3, 3))
    for j in range(3):
        e = np.zeros(3)
        e[j] = h
        J[:, :, j] = (func(pts + e) - func(pts - e)) / (2.0 * h)
    return J


def maxwell_checks(winding: Winding, sample_points: np.ndarray, h: Optional[float] = None,
                   ampere_loops: int = 3, seed: int = 0) -> Dict[str, float]:
    """Numerical verification of the vacuum Maxwell equations for a winding.

    * div_B_rel     max |∇·B| · L_char / |B|   (should be ~ finite-difference error)
    * curl_B_rel    max |∇×B| · L_char / |B| outside the conductor (should be ~0)
    * curlA_minus_B max |∇×A − B| / |B|
    * ampere_rel    |∮B·dl / (μ0 I) − 1| around loops encircling one conductor
    """

    pts = np.atleast_2d(np.asarray(sample_points, dtype=float))
    d = nearest_distance(pts, winding.points)
    pts = pts[d > 4.0 * winding.wire_radius_m]
    if h is None:
        h = 1e-3 * float(np.median(d))
    B = biot_savart(pts, winding)
    Bn = np.linalg.norm(B, axis=1)
    L = float(np.median(d))
    JB = _central_diff_field(lambda p: biot_savart(p, winding), pts, h)
    divB = JB[:, 0, 0] + JB[:, 1, 1] + JB[:, 2, 2]
    curlB = np.column_stack((JB[:, 2, 1] - JB[:, 1, 2], JB[:, 0, 2] - JB[:, 2, 0], JB[:, 1, 0] - JB[:, 0, 1]))
    JA = _central_diff_field(lambda p: vector_potential(p, winding), pts, h)
    curlA = np.column_stack((JA[:, 2, 1] - JA[:, 1, 2], JA[:, 0, 2] - JA[:, 2, 0], JA[:, 1, 0] - JA[:, 0, 1]))
    out = {
        "n_points": int(len(pts)),
        "fd_step_m": float(h),
        "div_B_rel_max": float(np.max(np.abs(divB) * L / Bn)),
        "curl_B_rel_max": float(np.max(np.linalg.norm(curlB, axis=1) * L / Bn)),
        "curlA_minus_B_rel_max": float(np.max(np.linalg.norm(curlA - B, axis=1) / Bn)),
    }
    # Ampère loops around the conductor at a few places along the active part
    rng = np.random.default_rng(seed)
    p1, dl = winding.segments()
    idx = rng.choice(np.arange(len(p1) // 8, len(p1) - len(p1) // 8), size=ampere_loops, replace=False)
    errs = []
    for k in idx:
        t = dl[k] / max(np.linalg.norm(dl[k]), 1e-300)
        trial = np.array([0, 0, 1.0]) if abs(t[2]) < 0.9 else np.array([1.0, 0, 0])
        n1 = trial - np.dot(trial, t) * t
        n1 /= np.linalg.norm(n1)
        n2 = np.cross(t, n1)
        rad = 3.0 * winding.wire_radius_m
        center = p1[k] + 0.5 * dl[k]
        th = np.linspace(0, 2 * np.pi, 721)
        loop = center + rad * (np.cos(th)[:, None] * n1 + np.sin(th)[:, None] * n2)
        Bl = biot_savart(loop, winding)
        tang = rad * (-np.sin(th)[:, None] * n1 + np.cos(th)[:, None] * n2)
        integrand = np.sum(Bl * tang, axis=1)
        circ = np.trapezoid(integrand, th)
        errs.append(abs(circ / (MU0 * winding.current_A)) - 1.0)
    out["ampere_rel_err_max"] = float(np.max(np.abs(errs)))
    out["ampere_loops"] = int(ampere_loops)
    return out


def return_sensitivity(winding: Winding, sample_points: np.ndarray) -> float:
    """RMS field change inside the vessel when the return path is dropped, relative to the full field."""

    full = biot_savart(sample_points, winding)
    active = Winding("active-only", winding.active_points, winding.current_A, winding.wire_radius_m)
    part = biot_savart(sample_points, active)
    return float(np.sqrt(np.mean(np.sum((full - part) ** 2, axis=1)) / np.mean(np.sum(full * full, axis=1))))


# ---------------------------------------------------------------------------
# Analytic references
# ---------------------------------------------------------------------------


def loop_field_analytic(points: np.ndarray, radius_m: float, current_A: float) -> np.ndarray:
    """Exact field of a circular loop (axis z, centre origin) via elliptic integrals."""

    pts = np.atleast_2d(points)
    x, y, z = pts[:, 0], pts[:, 1], pts[:, 2]
    rho = np.hypot(x, y)
    a = radius_m
    alpha2 = a * a + rho * rho + z * z - 2.0 * a * rho
    beta2 = a * a + rho * rho + z * z + 2.0 * a * rho
    beta = np.sqrt(beta2)
    k2 = 1.0 - alpha2 / beta2
    K = ellipk(k2)
    E = ellipe(k2)
    C = MU0 * current_A / math.pi
    with np.errstate(divide="ignore", invalid="ignore"):
        Brho = np.where(rho > 1e-12,
                        C * z / (2.0 * alpha2 * beta * rho) * ((a * a + rho * rho + z * z) * E - alpha2 * K), 0.0)
    Bz = C / (2.0 * alpha2 * beta) * ((a * a - rho * rho - z * z) * E + alpha2 * K)
    phi = np.arctan2(y, x)
    return np.column_stack((Brho * np.cos(phi), Brho * np.sin(phi), Bz))


def solenoid_axis_analytic(z: np.ndarray, radius_m: float, length_m: float, turns: int, current_A: float) -> np.ndarray:
    """On-axis field of a finite solenoid (uniform sheet current), tesla."""

    nI = turns * current_A / length_m
    zp = z + 0.5 * length_m
    zm = z - 0.5 * length_m
    return 0.5 * MU0 * nI * (zp / np.sqrt(zp * zp + radius_m**2) - zm / np.sqrt(zm * zm + radius_m**2))


def infinite_wire_analytic(rho: np.ndarray, current_A: float) -> np.ndarray:
    return MU0 * current_A / (2.0 * math.pi * rho)


def analytic_references(current_A: float = 1000.0, radius_m: float = 0.1) -> Dict[str, float]:
    """Compare segment Biot–Savart against closed forms.  Returns max relative errors."""

    from .geometry import circular_loop, solenoid, straight_wire

    out = {}
    # circular loop, off-axis points
    loop = Winding("loop", circular_loop(radius_m, n=2880), current_A, 1e-4)
    rng = np.random.default_rng(1)
    pts = rng.uniform(-2 * radius_m, 2 * radius_m, size=(400, 3))
    rho = np.hypot(pts[:, 0], pts[:, 1])
    pts = pts[np.abs(rho - radius_m) > 0.05 * radius_m]
    num = biot_savart(pts, loop)
    exact = loop_field_analytic(pts, radius_m, current_A)
    out["loop_offaxis_rel_max"] = float(np.max(np.linalg.norm(num - exact, axis=1) / np.linalg.norm(exact, axis=1)))
    # loop centre
    out["loop_centre_rel"] = float(abs(biot_savart(np.zeros((1, 3)), loop)[0, 2] / (MU0 * current_A / (2 * radius_m)) - 1.0))
    # finite solenoid on axis (helix vs ideal sheet: differs by the helical pitch term ~ 1/turns)
    turns, length = 200, 1.0
    sol = Winding("sol", solenoid(radius_m, length, turns, n_per_turn=120), current_A, 1e-4)
    z = np.linspace(-0.4 * length, 0.4 * length, 9)
    num_z = biot_savart(np.column_stack((np.zeros_like(z), np.zeros_like(z), z)), sol)[:, 2]
    out["solenoid_axis_rel_max"] = float(np.max(np.abs(num_z / solenoid_axis_analytic(z, radius_m, length, turns, current_A) - 1.0)))
    # long straight wire vs infinite-wire law at rho << L
    wire = Winding("wire", straight_wire(200.0, n=2), current_A, 1e-4)
    rho = np.array([0.01, 0.1, 1.0])
    num_w = biot_savart(np.column_stack((rho, np.zeros_like(rho), np.zeros_like(rho))), wire)[:, 1]
    out["straight_wire_rel_max"] = float(np.max(np.abs(num_w / infinite_wire_analytic(rho, current_A) - 1.0)))
    # interior of the conductor: |B| should rise linearly to the surface value
    thick = Winding("thick", straight_wire(200.0, n=2), current_A, 0.01)
    r_in = np.array([0.0025, 0.005, 0.0075, 0.01])
    b_in = np.abs(biot_savart(np.column_stack((r_in, np.zeros_like(r_in), np.zeros_like(r_in))), thick)[:, 1])
    expected = MU0 * current_A * r_in / (2 * math.pi * 0.01**2)
    out["conductor_interior_rel_max"] = float(np.max(np.abs(b_in / expected - 1.0)))
    return out
