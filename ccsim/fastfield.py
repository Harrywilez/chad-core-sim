"""Fast exact Biot–Savart (numba) for field-line tracing and guiding-centre pushes.

Same finite-segment Hanson–Hirshman kernel with uniform-current conductor interior as
`fields.biot_savart`; verified against it in exp01.  Falls back to the numpy version
when numba is not installed (slower, identical numbers).
"""

from __future__ import annotations

import math

import numpy as np

from .constants import MU0

try:
    from numba import njit, prange

    HAVE_NUMBA = True
except Exception:  # pragma: no cover
    HAVE_NUMBA = False

    def njit(*a, **k):
        def deco(f):
            return f
        return deco if not (a and callable(a[0])) else a[0]

    prange = range


class SegmentSource:
    """Flattened segment list of a Winding / MultiWinding: p1 (S,3), dl (S,3), I (S,), a (S,)."""

    def __init__(self, winding, current_scale: float = 1.0, coarsen: int = 1):
        parts = winding.windings if hasattr(winding, "windings") else [winding]
        P1, DL, I, A = [], [], [], []
        for w in parts:
            pts = w.points[::coarsen] if coarsen > 1 else w.points
            if coarsen > 1 and w.closed and not np.allclose(pts[-1], w.points[-1]):
                pts = np.vstack((pts, w.points[-1:]))
            p1 = pts[:-1]
            dl = pts[1:] - pts[:-1]
            L = np.linalg.norm(dl, axis=1)
            keep = L > 1e-15
            P1.append(p1[keep]); DL.append(dl[keep]); I.append(np.full(keep.sum(), w.current_A * current_scale)); A.append(np.full(keep.sum(), w.wire_radius_m))
        self.p1 = np.ascontiguousarray(np.vstack(P1))
        self.dl = np.ascontiguousarray(np.vstack(DL))
        self.I = np.ascontiguousarray(np.concatenate(I))
        self.a = np.ascontiguousarray(np.concatenate(A))
        self.L = np.linalg.norm(self.dl, axis=1)
        self.that = np.ascontiguousarray(self.dl / self.L[:, None])
        self.p2 = np.ascontiguousarray(self.p1 + self.dl)

    def B(self, points: np.ndarray) -> np.ndarray:
        pts = np.ascontiguousarray(np.atleast_2d(points), dtype=float)
        return _bs_kernel(pts, self.p1, self.p2, self.that, self.L, self.I, self.a)

    def B_and_grad(self, points: np.ndarray, h: float = 1e-4):
        """B and ∇B (M,3,3) by central differences: J[i, j] = ∂B_i/∂x_j."""
        pts = np.atleast_2d(points)
        M = len(pts)
        P = np.repeat(pts, 7, axis=0)
        for d in range(3):
            P[1 + 2 * d::7, d] += h
            P[2 + 2 * d::7, d] -= h
        B = self.B(P).reshape(M, 7, 3)
        J = np.empty((M, 3, 3))
        for d in range(3):
            J[:, :, d] = (B[:, 1 + 2 * d, :] - B[:, 2 + 2 * d, :]) / (2 * h)
        return B[:, 0, :], J


@njit(cache=True, parallel=False, fastmath=True)
def _bs_kernel(pts, p1, p2, that, L, I, a):
    M = pts.shape[0]
    S = p1.shape[0]
    out = np.zeros((M, 3))
    pref = MU0 / (4.0 * math.pi)
    for m in range(M):
        x, y, z = pts[m, 0], pts[m, 1], pts[m, 2]
        bx = 0.0; by = 0.0; bz = 0.0
        for s in range(S):
            Rix = x - p1[s, 0]; Riy = y - p1[s, 1]; Riz = z - p1[s, 2]
            Rfx = x - p2[s, 0]; Rfy = y - p2[s, 1]; Rfz = z - p2[s, 2]
            ri = math.sqrt(Rix * Rix + Riy * Riy + Riz * Riz)
            rf = math.sqrt(Rfx * Rfx + Rfy * Rfy + Rfz * Rfz)
            dot = Rix * Rfx + Riy * Rfy + Riz * Rfz
            rirf = ri * rf
            denom = rirf * (rirf + dot)
            if denom <= 1e-300:
                continue
            s0 = Rix * that[s, 0] + Riy * that[s, 1] + Riz * that[s, 2]
            rho2 = ri * ri - s0 * s0
            if rho2 < 0.0:
                rho2 = 0.0
            aw = a[s]
            factor = 1.0
            if rho2 < aw * aw and s0 >= -aw and s0 <= L[s] + aw:
                factor = rho2 / (aw * aw)
            coef = (ri + rf) / denom * factor * I[s]
            bx += (Riy * Rfz - Riz * Rfy) * coef
            by += (Riz * Rfx - Rix * Rfz) * coef
            bz += (Rix * Rfy - Riy * Rfx) * coef
        out[m, 0] = pref * bx; out[m, 1] = pref * by; out[m, 2] = pref * bz
    return out
