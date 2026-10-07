"""Rung 1: flux surfaces of a vacuum field — Poincaré sections, rotational transform,
magnetic well, island / chaos / open classification.

Field lines are traced on the *exact* Biot–Savart field (numba kernel) so that the
interpolation divergence error of a grid cannot fake or destroy surfaces.  All lines
are advanced together (batched RK4, fixed arclength step).

Definitions used
  * Poincaré section: punctures of the half-plane φ = 0, x > 0 (toroidal angle about z).
  * magnetic axis: centroid of the punctures of the line with the smallest puncture
    spread, refined by one pass of the return-map fixed-point search.
  * rotational transform ι = (poloidal angle advance about the axis) / (2π · toroidal turns).
  * surface test: punctures (r, θ) about the axis fitted with a Fourier series in θ
    (m ≤ 6); relative residual < `surface_tol` → nested surface; a line that stays in but
    fails the fit is 'island/chaotic'; a line that reaches the wall or a conductor is 'open'.
  * specific volume U = ∮ dl / B per toroidal turn (∝ V'(ψ)); well depth
    W = 1 − U_edge / U_axis  (> 0 is a magnetic well).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np
from scipy.spatial import cKDTree

from .fastfield import SegmentSource


@dataclass
class LineTrace:
    seed: np.ndarray
    punctures: np.ndarray          # (K, 2) = (R, z) at φ = 0 crossings
    turns: int                     # completed toroidal turns
    ended: str                     # 'complete' | 'wall' | 'wire' | 'stalled'
    length: float
    U_per_turn: np.ndarray         # ∮dl/B per toroidal turn
    mean_B: float
    kind: str = "open"
    iota: float = float("nan")
    r_mean: float = float("nan")
    residual: float = float("nan")


def torus_frame(active_points: np.ndarray):
    """Major radius and tube radius of the round torus the strands lie on (z-axis symmetric)."""
    rho = np.hypot(active_points[:, 0], active_points[:, 1])
    R0 = 0.5 * (rho.min() + rho.max())
    return float(R0), float(0.5 * (rho.max() - rho.min()))


def trace_batch(src: SegmentSource, seeds: np.ndarray, step: float, max_turns: int, wall: float,
                conductor_tree: Optional[cKDTree] = None, wire_clear: float = 0.008, max_steps: int = 60000,
                direction: float = 1.0) -> List[LineTrace]:
    x = np.array(seeds, dtype=float)
    n = len(x)
    alive = np.ones(n, bool)
    ended = np.array(["complete"] * n, dtype=object)
    phi_prev = np.arctan2(x[:, 1], x[:, 0])
    turns = np.zeros(n, int)
    punct = [[] for _ in range(n)]
    Uacc = np.zeros(n)
    Uturn = [[] for _ in range(n)]
    Bsum = np.zeros(n); Bcnt = np.zeros(n)
    length = np.zeros(n)

    def tangent(p):
        B = src.B(p)
        nb = np.linalg.norm(B, axis=1)
        nb = np.where(nb > 0, nb, 1e-300)
        return direction * B / nb[:, None], nb

    for it in range(max_steps):
        idx = np.where(alive)[0]
        if idx.size == 0:
            break
        p = x[idx]
        k1, nb = tangent(p)
        k2, _ = tangent(p + 0.5 * step * k1)
        k3, _ = tangent(p + 0.5 * step * k2)
        k4, _ = tangent(p + step * k3)
        pn = p + (step / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        Uacc[idx] += step / nb
        Bsum[idx] += nb; Bcnt[idx] += 1
        length[idx] += step
        phi = np.arctan2(pn[:, 1], pn[:, 0])
        # crossing of φ = 0 with φ increasing (going from negative to positive through 0), or decreasing for direction<0
        pp = phi_prev[idx]
        crossed = (pp < 0) & (phi >= 0) & (np.abs(phi - pp) < np.pi) if direction > 0 else (pp > 0) & (phi <= 0) & (np.abs(phi - pp) < np.pi)
        for j in np.where(crossed)[0]:
            i = idx[j]
            f = -pp[j] / (phi[j] - pp[j])
            q = p[j] + f * (pn[j] - p[j])
            punct[i].append((math.hypot(q[0], q[1]), q[2]))
            turns[i] += 1
            Uturn[i].append(Uacc[i]); Uacc[i] = 0.0
        phi_prev[idx] = phi
        x[idx] = pn
        r = np.linalg.norm(pn, axis=1)
        out = r >= wall
        if conductor_tree is not None and it % 5 == 0:
            d, _ = conductor_tree.query(pn)
            hit = d <= wire_clear
        else:
            hit = np.zeros(len(idx), bool)
        done = turns[idx] >= max_turns
        for j in range(len(idx)):
            if out[j]:
                ended[idx[j]] = "wall"; alive[idx[j]] = False
            elif hit[j]:
                ended[idx[j]] = "wire"; alive[idx[j]] = False
            elif done[j]:
                alive[idx[j]] = False
    for i in np.where(alive)[0]:
        ended[i] = "stalled"
    return [LineTrace(seeds[i], np.array(punct[i]).reshape(-1, 2), int(turns[i]), str(ended[i]), float(length[i]), np.array(Uturn[i]),
                      float(Bsum[i] / max(Bcnt[i], 1))) for i in range(n)]


def _fourier_residual(r: np.ndarray, th: np.ndarray, m_max: int = 6) -> float:
    cols = [np.ones_like(th)]
    for m in range(1, m_max + 1):
        cols += [np.cos(m * th), np.sin(m * th)]
    A = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(A, r, rcond=None)
    return float(np.sqrt(np.mean((A @ coef - r) ** 2)) / max(np.mean(r), 1e-12))


def classify(traces: List[LineTrace], axis_guess: np.ndarray, min_punctures: int = 12, surface_tol: float = 0.06):
    """Find the axis, compute ι and the surface test for every line.  Returns (axis, traces)."""
    # innermost surviving line → axis estimate
    cands = [t for t in traces if t.ended == "complete" and len(t.punctures) >= min_punctures]
    if cands:
        spread = [np.std(t.punctures, axis=0).sum() for t in cands]
        inner = cands[int(np.argmin(spread))]
        axis = inner.punctures.mean(axis=0)
    else:
        axis = np.asarray(axis_guess, dtype=float)
    for t in traces:
        P = t.punctures
        if len(P) < min_punctures or t.ended != "complete":
            t.kind = "open" if t.ended in ("wall", "wire") else "short"
            continue
        dR, dz = P[:, 0] - axis[0], P[:, 1] - axis[1]
        r = np.hypot(dR, dz)
        th = np.unwrap(np.arctan2(dz, dR))
        t.iota = float((th[-1] - th[0]) / (2 * np.pi * (len(P) - 1)))
        t.r_mean = float(r.mean())
        t.residual = _fourier_residual(r, np.mod(th, 2 * np.pi))
        t.kind = "surface" if t.residual < surface_tol else "island/chaotic"
    return axis, traces


def specific_volume(t: LineTrace) -> float:
    return float(np.mean(t.U_per_turn)) if len(t.U_per_turn) else float("nan")


def analyse_surfaces(winding, seeds: np.ndarray, axis_guess, step: float = 0.006, max_turns: int = 30, wall: float = 0.82,
                     wire_clear: float = 0.008, current_scale: float = 1.0) -> dict:
    src = SegmentSource(winding, current_scale)
    tree = cKDTree(winding.active_points[::2])
    traces = trace_batch(src, seeds, step, max_turns, wall, tree, wire_clear)
    axis, traces = classify(traces, axis_guess)
    surf = [t for t in traces if t.kind == "surface"]
    surf.sort(key=lambda t: t.r_mean)
    out = {
        "n_seeds": len(traces),
        "axis_R": float(axis[0]), "axis_z": float(axis[1]),
        "frac_surface": float(np.mean([t.kind == "surface" for t in traces])),
        "frac_island_chaotic": float(np.mean([t.kind == "island/chaotic" for t in traces])),
        "frac_open": float(np.mean([t.kind == "open" for t in traces])),
        "frac_short": float(np.mean([t.kind == "short" for t in traces])),
        "outermost_surface_r": float(surf[-1].r_mean) if surf else 0.0,
        "iota_axis": float(surf[0].iota) if surf else float("nan"),
        "iota_edge": float(surf[-1].iota) if surf else float("nan"),
        "iota_profile": [(float(t.r_mean), float(t.iota)) for t in surf],
        "well_depth": float(1.0 - specific_volume(surf[-1]) / specific_volume(surf[0])) if len(surf) >= 2 else float("nan"),
        "U_profile": [(float(t.r_mean), specific_volume(t)) for t in surf],
        "mean_B_axis": float(surf[0].mean_B) if surf else float("nan"),
        "per_line": [{"kind": t.kind, "iota": t.iota, "r": t.r_mean, "residual": t.residual, "ended": t.ended, "turns": t.turns} for t in traces],
        "punctures": [t.punctures.tolist() for t in traces],
    }
    return out


def seeds_on_ray(R0: float, r_tube: float, n: int, z0: float = 0.0, s_min: float = 0.05, s_max: float = 0.92) -> np.ndarray:
    s = np.linspace(s_min, s_max, n)
    return np.column_stack((R0 + s * r_tube, np.zeros(n), np.full(n, z0)))


def find_axis(winding, R0: float, r_tube: float, current_scale: float = 1.0, n_coarse: int = 16, turns: int = 8,
              wall: float = 1.05, wire_clear: float = 0.02) -> np.ndarray:
    """Magnetic axis at φ = 0: trace a coarse fan of seeds in the poloidal plane, take the
    surviving line with the smallest puncture spread, and use its puncture centroid.
    Refined once by re-seeding at that centroid."""

    src = SegmentSource(winding, current_scale)
    tree = cKDTree(winding.active_points[::2])
    rr = np.linspace(-0.55, 0.55, n_coarse) * r_tube
    seeds = np.column_stack((R0 + rr, np.zeros(n_coarse), np.zeros(n_coarse)))
    seeds = np.vstack((seeds, np.column_stack((np.full(n_coarse // 2, R0), np.zeros(n_coarse // 2), np.linspace(-0.4, 0.4, n_coarse // 2) * r_tube))))
    best = np.array([R0, 0.0])
    for _ in range(2):
        traces = trace_batch(src, seeds, 0.008, turns, wall, tree, wire_clear)
        cands = [t for t in traces if t.ended == "complete" and len(t.punctures) >= turns - 1]
        if not cands:
            break
        spread = [np.std(t.punctures, axis=0).sum() for t in cands]
        t = cands[int(np.argmin(spread))]
        best = t.punctures.mean(axis=0)
        seeds = np.column_stack((best[0] + np.linspace(-0.08, 0.08, 9) * r_tube, np.zeros(9), np.full(9, best[1])))
    return best


def analyse_from_axis(winding, R0: float, r_tube: float, n_seeds: int = 24, max_turns: int = 40, current_scale: float = 1.0,
                      wall: float = 1.05, wire_clear: float = 0.02, step: float = 0.006) -> dict:
    """Two-pass surface analysis: locate the magnetic axis, then seed a ray from it outboard
    up to the winding surface."""

    axis = find_axis(winding, R0, r_tube, current_scale, wall=wall, wire_clear=wire_clear)
    reach = max(R0 + r_tube - axis[0] - 0.03, 0.05)      # distance from the axis to the outboard strands
    s = np.linspace(0.03, 0.95, n_seeds)
    seeds = np.column_stack((axis[0] + s * reach, np.zeros(n_seeds), np.full(n_seeds, axis[1])))
    out = analyse_surfaces(winding, seeds, axis, step=step, max_turns=max_turns, wall=wall, wire_clear=wire_clear, current_scale=current_scale)
    out["axis_search"] = axis.tolist()
    out["reach"] = float(reach)
    return out
