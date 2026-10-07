"""Full-orbit charged-particle motion in SI units.

* Relativistic Boris pusher (reduces to the classical Boris rotation for
  v ≪ c); exact energy conservation in a pure magnetic field up to round-off.
* Per-particle adaptive sub-stepping keeps |Ω|·dt_sub ≤ `omega_dt_max`, so a
  particle that wanders into the strong field next to a conductor is still
  resolved while the bulk of the ensemble runs at the coarse step.
* Losses: spherical (or cylindrical) wall, conductor strike, leaving the field
  grid.  Every loss is time-stamped and classified.
* Diagnostics: kinetic energy, magnetic moment μ = m v⊥²/(2B), pitch, maximum
  radial excursion, adiabaticity ρ_L/L_B.

The field is supplied as an object with `B_at(x)` and `E_at(x, t)` methods;
`ScaledGridField` builds one from a unit-scale `FieldGrid` using the exact
Biot–Savart scaling B(s·x; I) = (I / I_grid) · (1/s) · B_grid(x) so a single
grid per geometry serves benchtop, reactor and dimensionless runs.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Callable, Dict, Optional

import numpy as np

from .constants import C_LIGHT, E_CHARGE, Species
from .fields import Electrostatics, FieldGrid
from .geometry import nearest_distance


# ---------------------------------------------------------------------------
# Field wrappers
# ---------------------------------------------------------------------------


class ScaledGridField:
    """Vacuum B from a unit-scale grid, plus optional induced and electrostatic E."""

    def __init__(self, grid: FieldGrid, scale: float, current_A: float,
                 electrostatics: Optional[Electrostatics] = None,
                 current_ramp: Optional[Callable[[float], float]] = None):
        self.grid = grid
        self.scale = float(scale)
        self.current_A = float(current_A)
        self.gain = self.current_A / (grid.winding.current_A * self.scale)
        self.electrostatics = electrostatics
        self.current_ramp = current_ramp          # g(t): I(t) = I * g(t); dI/dt from finite difference
        self.half_extent = grid.half * self.scale
        self.wire_radius = grid.winding.wire_radius_m * self.scale
        self._wire_points = grid.winding.points[::2] * self.scale

    def B_at(self, x: np.ndarray, t: float = 0.0) -> np.ndarray:
        g = self.gain * (self.current_ramp(t) if self.current_ramp else 1.0)
        return g * self.grid.B_at(x / self.scale)

    def E_at(self, x: np.ndarray, t: float = 0.0) -> np.ndarray:
        E = np.zeros_like(x)
        if self.current_ramp is not None:
            if self.grid.A is None:
                raise ValueError("induced E needs a grid built with with_A=True")
            h = 1e-9
            dgdt = (self.current_ramp(t + h) - self.current_ramp(t - h)) / (2 * h)
            # A scales as I (no 1/s): A(s x) = s * A_unit(x) * (I/I_grid) / s  -> A(s x; I) = (I/I_grid) A_grid(x)
            E += -(self.current_A / self.grid.winding.current_A) * dgdt * self.grid.A_at(x / self.scale)
        if self.electrostatics is not None:
            E += self.electrostatics.field(x)
        return E

    def Bmag_at(self, x: np.ndarray) -> np.ndarray:
        return abs(self.gain) * self.grid.Bmag_at(x / self.scale)

    def wire_distance_at(self, x: np.ndarray) -> np.ndarray:
        return self.scale * self.grid.wire_distance_at(x / self.scale)

    def exact_wire_distance(self, x: np.ndarray) -> np.ndarray:
        return nearest_distance(x, self._wire_points)

    def max_B(self) -> float:
        return float(abs(self.gain) * self.grid.Bmag.max())


# ---------------------------------------------------------------------------
# Particle ensemble
# ---------------------------------------------------------------------------


@dataclass
class Ensemble:
    species: Species
    x: np.ndarray                 # (N, 3) m
    v: np.ndarray                 # (N, 3) m/s
    alive: np.ndarray = None
    loss_time: np.ndarray = None
    loss_kind: np.ndarray = None
    weight: float = 1.0           # physical particles per marker (for rate ledgers)

    def __post_init__(self):
        n = len(self.x)
        # always copy: run_orbits integrates in place and callers routinely pass shared arrays
        self.x = np.array(self.x, dtype=float, copy=True)
        self.v = np.array(self.v, dtype=float, copy=True)
        if self.alive is None:
            self.alive = np.ones(n, dtype=bool)
        if self.loss_time is None:
            self.loss_time = np.full(n, np.inf)
        if self.loss_kind is None:
            self.loss_kind = np.array(["retained"] * n, dtype=object)

    @property
    def n(self) -> int:
        return len(self.x)

    def kinetic_energy_J(self) -> np.ndarray:
        v2 = np.sum(self.v * self.v, axis=1)
        gamma = 1.0 / np.sqrt(1.0 - np.minimum(v2 / C_LIGHT**2, 0.999999))
        # (gamma - 1) m c^2 written as m gamma^2 v^2 / (gamma + 1): exact and free of cancellation at v << c
        return self.species.mass_kg * gamma * gamma * v2 / (gamma + 1.0)

    def kinetic_energy_eV(self) -> np.ndarray:
        return self.kinetic_energy_J() / E_CHARGE


def maxwellian_ensemble(species: Species, n: int, temperature_eV: float, positions: np.ndarray, rng) -> Ensemble:
    sigma = math.sqrt(temperature_eV * E_CHARGE / species.mass_kg)
    v = rng.normal(0.0, sigma, size=(n, 3))
    return Ensemble(species, np.asarray(positions, dtype=float)[:n], v)


def monoenergetic_ensemble(species: Species, n: int, energy_eV: float, positions: np.ndarray, rng) -> Ensemble:
    speed = math.sqrt(2.0 * energy_eV * E_CHARGE / species.mass_kg)
    d = rng.normal(size=(n, 3))
    d /= np.linalg.norm(d, axis=1)[:, None]
    return Ensemble(species, np.asarray(positions, dtype=float)[:n], speed * d)


def sample_positions_in_vessel(n: int, r_min: float, r_max: float, rng, avoid_points: Optional[np.ndarray] = None,
                               min_distance: float = 0.0) -> np.ndarray:
    out = []
    while len(out) < n:
        c = rng.uniform(-r_max, r_max, size=(4 * n, 3))
        r = np.linalg.norm(c, axis=1)
        c = c[(r >= r_min) & (r <= r_max)]
        if avoid_points is not None and min_distance > 0:
            c = c[nearest_distance(c, avoid_points) >= min_distance]
        out.extend(c.tolist())
    return np.asarray(out[:n])


# ---------------------------------------------------------------------------
# Boris pusher
# ---------------------------------------------------------------------------


def _boris_step(x, v, q_over_m, B, E, dt):
    """Relativistic Boris (Boris–Vay-free classic form with γ at half steps)."""

    c2 = C_LIGHT**2
    # u = γ v
    v2 = np.sum(v * v, axis=1)
    gamma = 1.0 / np.sqrt(1.0 - np.minimum(v2 / c2, 0.999999))
    u = gamma[:, None] * v
    half_E = 0.5 * q_over_m * dt * E
    u_minus = u + half_E
    gamma_m = np.sqrt(1.0 + np.sum(u_minus * u_minus, axis=1) / c2)
    t = (0.5 * q_over_m * dt / gamma_m)[:, None] * B
    t2 = np.sum(t * t, axis=1)
    s = 2.0 * t / (1.0 + t2)[:, None]
    u_prime = u_minus + np.cross(u_minus, t)
    u_plus = u_minus + np.cross(u_prime, s)
    u_new = u_plus + half_E
    gamma_new = np.sqrt(1.0 + np.sum(u_new * u_new, axis=1) / c2)
    v_new = u_new / gamma_new[:, None]
    x_new = x + v_new * dt
    return x_new, v_new


@dataclass
class OrbitResult:
    t_final: float
    retained: int
    n: int
    retention: float
    wilson_lo: float
    wilson_hi: float
    counts: Dict[str, int]
    median_loss_time: float
    median_censored: bool
    mean_loss_time_restricted: float
    energy_drift_rel_max: float
    max_radius_p90: float
    mu_variation_median: float
    adiabaticity_median: float
    steps: int
    substep_max: int
    loss_positions: np.ndarray
    loss_kinds: np.ndarray
    loss_times: np.ndarray
    extra: dict = field(default_factory=dict)


def wilson(successes: int, n: int, z: float = 1.96):
    if n == 0:
        return 0.0, 0.0
    p = successes / n
    den = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4.0 * n * n)) / den
    return max(0.0, center - half), min(1.0, center + half)


def run_orbits(ens: Ensemble, field, t_max: float, dt: float, wall_radius: float,
               omega_dt_max: float = 0.2, wall_shape: str = "sphere", wall_half_length: Optional[float] = None,
               collision_operator: Optional[Callable] = None, event_hook: Optional[Callable] = None,
               sample_every: int = 0, max_substeps: int = 4096) -> OrbitResult:
    """Integrate an ensemble until t_max, recording losses.

    `collision_operator(ens, alive_ids, dt, t)` may modify velocities in place
    (Coulomb / neutral collisions).  `event_hook(ens, alive_ids, dt, t)` is
    called each coarse step (e.g. fusion event sampling) and may return a dict
    of counters that are summed into `extra`.
    """

    q_over_m = ens.species.charge_C / ens.species.mass_kg
    x, v = ens.x, ens.v
    n = ens.n
    E0 = ens.kinetic_energy_J().copy()
    r0 = np.linalg.norm(x, axis=1)
    max_r = r0.copy()
    B_init = field.B_at(x)
    Bn0 = np.linalg.norm(B_init, axis=1)
    vpar0 = np.sum(v * B_init, axis=1) / np.maximum(Bn0, 1e-300)
    vperp2_0 = np.maximum(np.sum(v * v, axis=1) - vpar0**2, 0.0)
    mu0 = ens.species.mass_kg * vperp2_0 / (2.0 * np.maximum(Bn0, 1e-300))
    mu_min = mu0.copy()
    mu_max = mu0.copy()
    energy_drift = 0.0
    substep_max = 1
    steps = int(math.ceil(t_max / dt))
    t = 0.0
    extra: Dict[str, float] = {}
    samples = []
    has_E = getattr(field, "current_ramp", None) is not None or getattr(field, "electrostatics", None) is not None
    track_energy = collision_operator is None and not has_E
    for step in range(steps):
        ids = np.flatnonzero(ens.alive)
        if len(ids) == 0:
            break
        xi, vi = x[ids], v[ids]
        B = field.B_at(xi, t)
        Bn = np.linalg.norm(B, axis=1)
        omega_dt = abs(q_over_m) * Bn * dt
        nsub = np.minimum(np.maximum(np.ceil(omega_dt / omega_dt_max).astype(int), 1), max_substeps)
        substep_max = max(substep_max, int(nsub.max()))
        for k in np.unique(nsub):
            sel = np.flatnonzero(nsub == k)
            xs, vs = xi[sel], vi[sel]
            sub_dt = dt / k
            if k == 1:
                E = field.E_at(xs, t)
                xs, vs = _boris_step(xs, vs, q_over_m, B[sel], E, sub_dt)
            else:
                tt = t
                for _ in range(int(k)):
                    Bs = field.B_at(xs, tt)
                    Es = field.E_at(xs, tt)
                    xs, vs = _boris_step(xs, vs, q_over_m, Bs, Es, sub_dt)
                    tt += sub_dt
            xi[sel], vi[sel] = xs, vs
        x[ids], v[ids] = xi, vi
        t += dt
        if collision_operator is not None:
            collision_operator(ens, ids, dt, t)
        if event_hook is not None:
            counters = event_hook(ens, ids, dt, t)
            if counters:
                for key, val in counters.items():
                    extra[key] = extra.get(key, 0.0) + val
        # diagnostics
        xi, vi = x[ids], v[ids]
        r = np.linalg.norm(xi, axis=1)
        max_r[ids] = np.maximum(max_r[ids], r)
        if track_energy:
            Ek = ens.kinetic_energy_J()[ids]
            energy_drift = max(energy_drift, float(np.max(np.abs(Ek - E0[ids]) / E0[ids])))
        Bnow = field.B_at(xi, t)
        Bn = np.linalg.norm(Bnow, axis=1)
        vpar = np.sum(vi * Bnow, axis=1) / np.maximum(Bn, 1e-300)
        vperp2 = np.maximum(np.sum(vi * vi, axis=1) - vpar**2, 0.0)
        mu = ens.species.mass_kg * vperp2 / (2.0 * np.maximum(Bn, 1e-300))
        mu_min[ids] = np.minimum(mu_min[ids], mu)
        mu_max[ids] = np.maximum(mu_max[ids], mu)
        # losses
        if wall_shape == "sphere":
            wall = r >= wall_radius
        else:
            wall = (np.hypot(xi[:, 0], xi[:, 1]) >= wall_radius) | (np.abs(xi[:, 2]) >= (wall_half_length or wall_radius))
        outside = ~np.all(np.abs(xi) < field.half_extent, axis=1)
        approx = field.wire_distance_at(np.clip(xi, -field.half_extent * 0.999, field.half_extent * 0.999))
        maybe = approx <= max(3.0 * field.grid.dx * field.scale, 2.0 * field.wire_radius)
        wire = np.zeros(len(ids), dtype=bool)
        if np.any(maybe):
            wire[maybe] = field.exact_wire_distance(xi[maybe]) <= field.wire_radius
        bad = ~np.isfinite(xi).all(axis=1)
        lost = wall | outside | wire | bad
        if np.any(lost):
            li = ids[lost]
            ens.alive[li] = False
            ens.loss_time[li] = t
            kinds = np.where(bad[lost], "numerical", np.where(wire[lost], "wire", np.where(wall[lost], "wall", "grid-exit")))
            ens.loss_kind[li] = kinds
        if sample_every and (step + 1) % sample_every == 0:
            samples.append((t, x.copy(), ens.alive.copy()))
    retained = int(np.sum(ens.alive))
    lo, hi = wilson(retained, n)
    lt = np.where(np.isfinite(ens.loss_time), ens.loss_time, t)
    counts = {k: int(np.sum(ens.loss_kind == k)) for k in ("retained", "wall", "wire", "grid-exit", "numerical")}
    censored = retained > n / 2
    median = float(t) if censored else float(np.median(lt))
    mu_var = np.where(mu_max > 0, (mu_max - mu_min) / np.maximum(mu_max, 1e-300), 0.0)
    # adiabaticity: gyroradius over field-gradient scale, estimated at the initial positions
    h = 1e-3 * field.half_extent
    grad = np.zeros(n)
    for j in range(3):
        e = np.zeros(3)
        e[j] = h
        grad += ((field.Bmag_at(ens.x + e) - field.Bmag_at(ens.x - e)) / (2 * h)) ** 2
    L_B = Bn0 / np.maximum(np.sqrt(grad), 1e-300)
    rho_L = ens.species.mass_kg * np.sqrt(vperp2_0) / (abs(ens.species.charge_C) * np.maximum(Bn0, 1e-300))
    lost_mask = ~ens.alive
    return OrbitResult(
        t_final=float(t), retained=retained, n=n, retention=retained / n, wilson_lo=lo, wilson_hi=hi, counts=counts,
        median_loss_time=median, median_censored=censored, mean_loss_time_restricted=float(np.mean(lt)),
        energy_drift_rel_max=float(energy_drift), max_radius_p90=float(np.percentile(max_r, 90)),
        mu_variation_median=float(np.median(mu_var)), adiabaticity_median=float(np.median(rho_L / L_B)),
        steps=steps, substep_max=substep_max,
        loss_positions=ens.x[lost_mask].copy(), loss_kinds=ens.loss_kind[lost_mask].copy(), loss_times=ens.loss_time[lost_mask].copy(),
        extra={"samples": samples, **extra},
    )
