"""Exact Eckart-barrier transmission and thermal tunnelling factors.

A one-dimensional barrier is specified by the three numbers a converged
CI-NEB plus a saddle-point frequency calculation returns:
    V0   barrier height (eV, relative to reactants)
    dE   reaction energy (eV, products minus reactants; must be < V0)
    w    imaginary frequency magnitude at the saddle (cm^-1)
The Eckart form V(y) = dE·y/(1+y) + B·y/(1+y)², y = e^{2x/L}, matches V0, dE
and the saddle curvature; its transmission T(E) is closed-form (Eckart 1930).
Γ(T) = ⟨T(E)⟩_Boltzmann / e^{−V0/kT} multiplies an Arrhenius rate.
"""

from __future__ import annotations

import math

import numpy as np

from .constants import HBAR, K_B, E_CHARGE

HC_EV_CM = 1.2398419843320026e-4


def _logcosh(x):
    ax = np.abs(np.asarray(x, dtype=float))
    return ax + np.log1p(np.exp(-2.0 * ax)) - math.log(2.0)


def eckart_transmission(E_eV, V0_eV: float, dE_eV: float, imag_freq_cm1: float) -> np.ndarray:
    """Exact transmission probability through an asymmetric Eckart barrier."""

    E = np.atleast_1d(np.asarray(E_eV, dtype=float))
    hw = imag_freq_cm1 * HC_EV_CM
    A = dE_eV
    B = 2.0 * V0_eV - dE_eV + 2.0 * math.sqrt(V0_eV * (V0_eV - dE_eV))
    # Curvature: |V''| = 2B (2 s*(1-s*)/w)^2 with s* = (B+dE)/(2B)  ->  C = hbar^2/(2 m w^2) follows from hw
    s_star = (B + A) / (2.0 * B)
    kappa = 2.0 * B * (2.0 * s_star * (1.0 - s_star)) ** 2       # = |V''| * w^2  (eV)
    # hw^2 = 2 * C * kappa  (since hbar^2/(2m) V''/... )  -> C = hw^2 / (2 kappa)
    C = hw * hw / (2.0 * kappa)
    out = np.zeros_like(E)
    ok = (E > 0.0) & (E > A)
    e = E[ok]
    a = 0.5 * np.sqrt(e / C)
    b = 0.5 * np.sqrt((e - A) / C)
    d_term = math.cosh(math.pi * math.sqrt((B - C) / C)) if B >= C else math.cos(math.pi * math.sqrt((C - B) / C))
    log_d = math.log(d_term) if d_term > 0 else -math.inf
    num = np.logaddexp(_logcosh(2 * math.pi * (a - b)), log_d)
    den = np.logaddexp(_logcosh(2 * math.pi * (a + b)), log_d)
    out[ok] = 1.0 - np.exp(num - den)
    return np.clip(out, 0.0, 1.0)


def crossover_temperature_K(imag_freq_cm1: float) -> float:
    return imag_freq_cm1 * HC_EV_CM * E_CHARGE / (2.0 * math.pi * K_B)


def tunnelling_factor(T_K: float, V0_eV: float, dE_eV: float, imag_freq_cm1: float, n: int = 40000) -> float:
    """Γ(T) = ∫ e^{-x} T(x kT) dx / e^{-V0/kT}  (x = E / kT)."""

    kT = K_B * T_K / E_CHARGE
    x_max = V0_eV / kT + 60.0
    x = np.concatenate([np.linspace(0.0, min(5.0, x_max), 4001), np.linspace(min(5.0, x_max), x_max, n)])
    x = np.unique(x)
    t = eckart_transmission(x * kT, V0_eV, dE_eV, imag_freq_cm1)
    avg = float(np.trapezoid(np.exp(-x) * t, x))
    return avg / math.exp(-V0_eV / kT)


def imag_freq_from_width(V0_eV: float, dE_eV: float, width_A: float, mass_amu: float) -> float:
    """Q-Surface convention: width w (Å) and effective mass m (amu) → ω‡ (cm^-1); only m·w² matters."""

    HBAR2_OVER_2AMU_EV_A2 = 0.0020900796402483607
    B = 2.0 * V0_eV - dE_eV + 2.0 * math.sqrt(V0_eV * (V0_eV - dE_eV))
    s_star = (B + dE_eV) / (2.0 * B)
    curvature = 2.0 * B * (2.0 * s_star * (1.0 - s_star) / width_A) ** 2
    hw = math.sqrt(2.0 * HBAR2_OVER_2AMU_EV_A2 * curvature / mass_amu)
    return hw / HC_EV_CM
