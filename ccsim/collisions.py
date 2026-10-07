"""Test-particle Coulomb collisions against Maxwellian background species.

Implements the NRL Plasma Formulary relaxation rates for a test particle α
moving through field species β (slowing-down ν_s, transverse diffusion ν_⊥,
parallel diffusion ν_∥) as a Monte Carlo Langevin operator applied every
coarse time step.  This is the standard test-particle Fokker–Planck treatment:
the background is fixed (no self-consistency), collisions are binary and
small-angle (Coulomb logarithm ≫ 1), and the operator is exact in the limit
ν·dt ≪ 1, which `apply` enforces by sub-stepping.

Neutral (charge-exchange / elastic) collisions are represented by a constant
cross-section hard-sphere model — a placeholder to be replaced by tabulated
data (LXCat / Phelps) if neutral-dominated regimes matter.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import List, Optional, Sequence

import numpy as np
from scipy.special import erf

from .constants import E_CHARGE, EPS0, K_B, Species


@dataclass(frozen=True)
class Background:
    species: Species
    density_m3: float
    temperature_eV: float

    @property
    def thermal_speed(self) -> float:      # sqrt(2 k T / m)
        return math.sqrt(2.0 * self.temperature_eV * E_CHARGE / self.species.mass_kg)


def coulomb_log_ion_ion(test: Species, bg: Background, test_energy_eV: float) -> float:
    """NRL ion–ion Coulomb logarithm (n in cm^-3, T in eV); test particle at its energy-equivalent T."""

    n_cm3 = bg.density_m3 * 1e-6
    z1, z2 = abs(test.charge_number), abs(bg.species.charge_number)
    mu1, mu2 = test.mass_amu, bg.species.mass_amu
    T1 = max(test_energy_eV * 2.0 / 3.0, 1e-3)
    T2 = max(bg.temperature_eV, 1e-3)
    arg = (z1 * z2 * (mu1 + mu2) / (mu1 * T2 + mu2 * T1)) * math.sqrt(n_cm3 * z2 * z2 / T2)
    return max(2.0, 23.0 - math.log(max(arg, 1e-300)))


def coulomb_log_electron(bg_e: Background) -> float:
    """NRL electron–ion Coulomb logarithm for T_e > 10 Z² eV (n in cm^-3)."""

    n_cm3 = bg_e.density_m3 * 1e-6
    Te = max(bg_e.temperature_eV, 1e-3)
    return max(2.0, 24.0 - math.log(math.sqrt(n_cm3) / Te))


def _psi(x: np.ndarray) -> np.ndarray:
    # psi(x) = (2/sqrt(pi)) int_0^x sqrt(t) e^{-t} dt = erf(sqrt x) - (2/sqrt pi) sqrt(x) e^{-x}
    sx = np.sqrt(x)
    return erf(sx) - (2.0 / math.sqrt(math.pi)) * sx * np.exp(-x)


def _psi_prime(x: np.ndarray) -> np.ndarray:
    return (2.0 / math.sqrt(math.pi)) * np.sqrt(x) * np.exp(-x)


def relaxation_rates(test: Species, v: np.ndarray, bg: Background, ln_lambda: Optional[float] = None):
    """Return (nu_s, nu_perp, nu_par) [1/s] for test-particle speeds v (m/s).  NRL formulary."""

    speed = np.maximum(np.linalg.norm(v, axis=1), 1.0)
    if ln_lambda is None:
        E_eV = 0.5 * test.mass_kg * np.mean(speed) ** 2 / E_CHARGE
        ln_lambda = coulomb_log_electron(bg) if bg.species.name.startswith("e") else coulomb_log_ion_ion(test, bg, E_eV)
    nu0 = (bg.density_m3 * test.charge_C**2 * bg.species.charge_C**2 * ln_lambda
           / (4.0 * math.pi * EPS0**2 * test.mass_kg**2 * speed**3))
    x = bg.species.mass_kg * speed**2 / (2.0 * bg.temperature_eV * E_CHARGE)
    psi = _psi(x)
    psip = _psi_prime(x)
    nu_s = (1.0 + test.mass_kg / bg.species.mass_kg) * psi * nu0
    nu_perp = 2.0 * ((1.0 - 1.0 / (2.0 * x)) * psi + psip) * nu0
    nu_par = (psi / x) * nu0
    return nu_s, nu_perp, nu_par


def slowing_down_time(test: Species, energy_eV: float, bg: Background) -> float:
    v = np.array([[math.sqrt(2 * energy_eV * E_CHARGE / test.mass_kg), 0.0, 0.0]])
    nu_s, _, _ = relaxation_rates(test, v, bg)
    return float(1.0 / nu_s[0])


class CoulombOperator:
    """Langevin Monte Carlo operator: apply(ens, ids, dt, t) modifies velocities in place."""

    def __init__(self, backgrounds: Sequence[Background], rng, max_nu_dt: float = 0.1):
        self.backgrounds = list(backgrounds)
        self.rng = rng
        self.max_nu_dt = max_nu_dt

    def __call__(self, ens, ids, dt, t):
        v = ens.v[ids]
        for bg in self.backgrounds:
            nu_s, nu_perp, nu_par = relaxation_rates(ens.species, v, bg)
            nu_max = float(np.max(np.maximum(nu_s, np.maximum(nu_perp, nu_par)))) if len(v) else 0.0
            nsub = max(1, int(math.ceil(nu_max * dt / self.max_nu_dt)))
            sub = dt / nsub
            for _ in range(nsub):
                if nsub > 1:
                    nu_s, nu_perp, nu_par = relaxation_rates(ens.species, v, bg)
                speed = np.linalg.norm(v, axis=1)
                e1 = v / np.maximum(speed, 1e-300)[:, None]
                trial = np.where(np.abs(e1[:, 2:3]) < 0.9, np.array([[0.0, 0.0, 1.0]]), np.array([[1.0, 0.0, 0.0]]))
                e2 = trial - np.sum(trial * e1, axis=1)[:, None] * e1
                e2 /= np.maximum(np.linalg.norm(e2, axis=1), 1e-300)[:, None]
                e3 = np.cross(e1, e2)
                d_par = -nu_s * speed * sub + np.sqrt(np.maximum(nu_par, 0.0) * speed**2 * sub) * self.rng.normal(size=len(v))
                sig_perp = np.sqrt(np.maximum(nu_perp, 0.0) * speed**2 * sub / 2.0)
                d2 = sig_perp * self.rng.normal(size=len(v))
                d3 = sig_perp * self.rng.normal(size=len(v))
                v = v + d_par[:, None] * e1 + d2[:, None] * e2 + d3[:, None] * e3
        ens.v[ids] = v


class NeutralCollisionOperator:
    """Hard-sphere collisions with a neutral gas (placeholder cross-section).

    Each step, a particle collides with probability 1 − exp(−n σ v dt); on a
    collision its velocity is redrawn isotropically at the same speed (elastic,
    heavy target) or, with probability `charge_exchange_fraction`, replaced by a
    thermal neutral velocity (charge exchange).  Replace σ with tabulated data
    for quantitative work.
    """

    def __init__(self, density_m3: float, temperature_eV: float, target_mass_kg: float,
                 cross_section_m2: float = 3e-19, charge_exchange_fraction: float = 0.5, rng=None):
        self.n = density_m3
        self.T = temperature_eV
        self.m_t = target_mass_kg
        self.sigma = cross_section_m2
        self.cx = charge_exchange_fraction
        self.rng = rng or np.random.default_rng(0)
        self.count = 0

    def mean_free_path(self) -> float:
        return 1.0 / (self.n * self.sigma)

    def __call__(self, ens, ids, dt, t):
        v = ens.v[ids]
        speed = np.linalg.norm(v, axis=1)
        p = 1.0 - np.exp(-self.n * self.sigma * speed * dt)
        hit = self.rng.random(len(ids)) < p
        if not np.any(hit):
            return
        k = int(np.sum(hit))
        self.count += k
        d = self.rng.normal(size=(k, 3))
        d /= np.linalg.norm(d, axis=1)[:, None]
        cx = self.rng.random(k) < self.cx
        sigma_t = math.sqrt(self.T * E_CHARGE / ens.species.mass_kg)
        new = np.where(cx[:, None], self.rng.normal(0.0, sigma_t, size=(k, 3)), speed[hit][:, None] * d)
        v[hit] = new
        ens.v[ids] = v
