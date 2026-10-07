"""Rung 2: guiding-centre orbits with Coulomb pitch-angle scattering on the exact vacuum field.

Non-relativistic guiding-centre equations in a vacuum field (curvature κ = ∇⊥B/B):
    dX/dt  = v∥ b + (v∥² + v⊥²/2) / Ω · (b × ∇B) / B
    dv∥/dt = −(μ/m) b·∇B,        μ = m v⊥² / (2B)   (conserved)
Ω = qB/m.  Time-stepped with RK4 on (X, v∥); B and ∇B from the numba Biot–Savart
kernel by central differences.  Collisions: Lorentz pitch-angle scattering with rate ν
(Δξ = −ν ξ dt + √((1−ξ²) ν dt) N(0,1), ξ = v∥/v), energy fixed — the operator that
matters for trapped/passing transitions; slowing-down is left out on purpose (it only
rescales v).  ν can be scaled up to fit collision times into a tractable run; the
1/ν-regime extrapolation is done by the experiment script, not here.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .fastfield import SegmentSource


@dataclass
class GCResult:
    t: np.ndarray                # sample times
    alive_frac: np.ndarray       # survival curve
    loss_times: np.ndarray       # per particle (inf if retained)
    r_final: np.ndarray          # final minor-radius-like coordinate
    n_steps: int
    energy_drift_rel_max: float


def gc_rhs(src: SegmentSource, X: np.ndarray, vpar: np.ndarray, mu: np.ndarray, q_over_m: float, mass: float, h: float):
    B, J = src.B_and_grad(X, h)
    Bm = np.linalg.norm(B, axis=1)
    b = B / Bm[:, None]
    gradB = np.einsum("nij,ni->nj", J, b)          # ∇|B| = (J^T b)
    Omega = q_over_m * Bm
    vperp2 = 2.0 * mu * Bm / mass
    vd = ((vpar**2 + 0.5 * vperp2) / (Omega * Bm))[:, None] * np.cross(b, gradB)
    dX = vpar[:, None] * b + vd
    dvpar = -(mu / mass) * np.einsum("ni,ni->n", b, gradB)
    return dX, dvpar, Bm


def run_gc(src: SegmentSource, X0: np.ndarray, v0: np.ndarray, species, t_max: float, dt: float, loss_fn, nu: float = 0.0,
           rng=None, sample_every: int = 20, fd_h: float = 1e-4) -> GCResult:
    """X0 (N,3) positions (metres), v0 (N,3) velocities; loss_fn(X) -> bool array (lost)."""

    q_over_m = species.charge_C / species.mass_kg
    m = species.mass_kg
    N = len(X0)
    X = np.array(X0, dtype=float)
    B0 = src.B(X)
    Bm0 = np.linalg.norm(B0, axis=1)
    b0 = B0 / Bm0[:, None]
    vpar = np.einsum("ni,ni->n", v0, b0)
    vperp2 = np.sum(v0**2, axis=1) - vpar**2
    mu = m * vperp2 / (2.0 * Bm0)
    E0 = 0.5 * m * np.sum(v0**2, axis=1)
    alive = np.ones(N, bool)
    loss_times = np.full(N, np.inf)
    rng = rng or np.random.default_rng(0)
    n_steps = int(math.ceil(t_max / dt))
    ts, af = [], []
    e_drift = 0.0
    for step in range(n_steps):
        idx = np.where(alive)[0]
        if idx.size == 0:
            break
        x, vp, mu_i = X[idx], vpar[idx], mu[idx]
        k1x, k1v, _ = gc_rhs(src, x, vp, mu_i, q_over_m, m, fd_h)
        k2x, k2v, _ = gc_rhs(src, x + 0.5 * dt * k1x, vp + 0.5 * dt * k1v, mu_i, q_over_m, m, fd_h)
        k3x, k3v, _ = gc_rhs(src, x + 0.5 * dt * k2x, vp + 0.5 * dt * k2v, mu_i, q_over_m, m, fd_h)
        k4x, k4v, Bm = gc_rhs(src, x + dt * k3x, vp + dt * k3v, mu_i, q_over_m, m, fd_h)
        X[idx] = x + (dt / 6.0) * (k1x + 2 * k2x + 2 * k3x + k4x)
        vpar[idx] = vp + (dt / 6.0) * (k1v + 2 * k2v + 2 * k3v + k4v)
        if nu > 0:
            # Lorentz pitch-angle scattering at fixed speed
            B = src.B(X[idx]); Bm = np.linalg.norm(B, axis=1)
            v2 = vpar[idx] ** 2 + 2.0 * mu[idx] * Bm / m
            v = np.sqrt(v2)
            xi = np.clip(vpar[idx] / v, -1, 1)
            xi = xi - nu * xi * dt + np.sqrt(np.maximum(1 - xi**2, 0) * nu * dt) * rng.normal(size=idx.size)
            xi = np.clip(xi, -1, 1)
            vpar[idx] = xi * v
            mu[idx] = m * v2 * (1 - xi**2) / (2.0 * Bm)
        lost = loss_fn(X[idx])
        if lost.any():
            li = idx[lost]
            alive[li] = False
            loss_times[li] = (step + 1) * dt
        if (step + 1) % sample_every == 0:
            ts.append((step + 1) * dt); af.append(alive.mean())
            if idx.size and step % (sample_every * 10) == 0:
                B = src.B(X[idx]); Bm = np.linalg.norm(B, axis=1)
                E = 0.5 * m * (vpar[idx] ** 2 + 2.0 * mu[idx] * Bm / m)
                e_drift = max(e_drift, float(np.max(np.abs(E / E0[idx] - 1))))
    return GCResult(np.array(ts), np.array(af), loss_times, np.linalg.norm(X, axis=1), n_steps, e_drift)
