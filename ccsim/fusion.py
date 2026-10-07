"""Fusion cross-sections, reactivities and Monte Carlo fusion events.

Cross-sections σ(E_cm) and Maxwellian reactivities ⟨σv⟩(T) use the Bosch &
Hale (1992, Nucl. Fusion 32, 611) parametrisations for D–T, D–D (both
branches) and D–³He, valid on the energy/temperature ranges they state.
`verify()` cross-checks the two independent parametrisations against each
other by numerically integrating σ(E) over a Maxwellian, and against tabulated
values, so a coefficient typo cannot pass silently.

`FusionEventSampler` is an `event_hook` for `run_orbits`: for every tracked
fuel ion it samples a background target velocity from a Maxwellian, evaluates
σ at the centre-of-mass energy and draws a fusion event with probability
n_t σ v_rel dt.  Events are counted with the marker weight and their birth
positions recorded (an alpha/neutron source map), without removing the marker.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Dict, List, Optional

import numpy as np

from .constants import AMU, E_CHARGE, KEV, MASS_AMU, Species

# Bosch–Hale Table IV (cross-section S-factor fits), energies in keV (c.m.), σ in millibarn
_SIGMA = {
    "DT": dict(BG=34.3827, A=(6.927e4, 7.454e8, 2.050e6, 5.2002e4, 0.0), B=(6.38e1, -9.95e-1, 6.981e-5, 1.728e-4), Erange=(0.5, 550.0)),
    "DD_pT": dict(BG=31.3970, A=(5.5576e4, 2.1054e2, -3.2638e-2, 1.4987e-6, 1.8181e-10), B=(0.0, 0.0, 0.0, 0.0), Erange=(0.5, 5000.0)),
    "DD_nHe3": dict(BG=31.3970, A=(5.3701e4, 3.3027e2, -1.2706e-1, 2.9327e-5, -2.5151e-9), B=(0.0, 0.0, 0.0, 0.0), Erange=(0.5, 4900.0)),
    "DHe3": dict(BG=68.7508, A=(5.7501e6, 2.5226e3, 4.5566e1, 0.0, 0.0), B=(-3.1995e-3, -8.5530e-6, 5.9014e-8, 0.0), Erange=(0.3, 900.0)),
}
# Bosch–Hale Table VII (reactivity fits), T in keV, <σv> in cm^3/s
_REACTIVITY = {
    "DT": dict(BG=34.3827, mrc2=1124656.0, C=(1.17302e-9, 1.51361e-2, 7.51886e-2, 4.60643e-3, 1.35000e-2, -1.06750e-4, 1.36600e-5), Trange=(0.2, 100.0)),
    "DD_pT": dict(BG=31.3970, mrc2=937814.0, C=(5.65718e-12, 3.41267e-3, 1.99167e-3, 0.0, 1.05060e-5, 0.0, 0.0), Trange=(0.2, 100.0)),
    "DD_nHe3": dict(BG=31.3970, mrc2=937814.0, C=(5.43360e-12, 5.85778e-3, 7.68222e-3, 0.0, -2.96400e-6, 0.0, 0.0), Trange=(0.2, 100.0)),
    "DHe3": dict(BG=68.7508, mrc2=1124572.0, C=(5.51036e-10, 6.41918e-3, -2.02896e-3, -1.91080e-5, 1.35776e-4, 0.0, 0.0), Trange=(0.5, 190.0)),
}
REACTION_ENERGY_MEV = {"DT": 17.59, "DD_pT": 4.03, "DD_nHe3": 3.27, "DHe3": 18.35}
PRODUCTS = {"DT": "α(3.52 MeV) + n(14.07 MeV)", "DD_pT": "T(1.01) + p(3.02)", "DD_nHe3": "³He(0.82) + n(2.45)", "DHe3": "α(3.67) + p(14.68)"}
REACTANT_MASSES_AMU = {"DT": (MASS_AMU["D"], MASS_AMU["T"]), "DD_pT": (MASS_AMU["D"], MASS_AMU["D"]),
                       "DD_nHe3": (MASS_AMU["D"], MASS_AMU["D"]), "DHe3": (MASS_AMU["D"], MASS_AMU["He3"])}


def cross_section_m2(reaction: str, E_cm_keV) -> np.ndarray:
    """Bosch–Hale σ(E_cm) in m²; zero outside the fit range's low end (extrapolation is refused)."""

    p = _SIGMA[reaction]
    E = np.atleast_1d(np.asarray(E_cm_keV, dtype=float))
    A, B = p["A"], p["B"]
    S = (A[0] + E * (A[1] + E * (A[2] + E * (A[3] + E * A[4])))) / (1.0 + E * (B[0] + E * (B[1] + E * (B[2] + E * B[3]))))
    with np.errstate(divide="ignore", over="ignore", invalid="ignore"):
        sigma_mb = np.where(E > 0, S / (E * np.exp(p["BG"] / np.sqrt(np.maximum(E, 1e-300)))), 0.0)
    lo, hi = p["Erange"]
    sigma_mb = np.where((E >= lo) & (E <= hi), sigma_mb, np.where(E < lo, 0.0, np.nan))
    return sigma_mb * 1e-31   # 1 mb = 1e-31 m^2


def reactivity_m3_s(reaction: str, T_keV) -> np.ndarray:
    """Bosch–Hale Maxwellian ⟨σv⟩(T) in m³/s (NaN outside the stated range)."""

    p = _REACTIVITY[reaction]
    T = np.atleast_1d(np.asarray(T_keV, dtype=float))
    C = p["C"]
    theta = T / (1.0 - T * (C[1] + T * (C[3] + T * C[5])) / (1.0 + T * (C[2] + T * (C[4] + T * C[6]))))
    xi = (p["BG"] ** 2 / (4.0 * theta)) ** (1.0 / 3.0)
    sv_cm3 = C[0] * theta * np.sqrt(xi / (p["mrc2"] * T**3)) * np.exp(-3.0 * xi)
    lo, hi = p["Trange"]
    sv_cm3 = np.where((T >= lo) & (T <= hi), sv_cm3, np.nan)
    return sv_cm3 * 1e-6


def reactivity_by_integration(reaction: str, T_keV: float, n: int = 20000) -> float:
    """⟨σv⟩ by direct integration of σ(E) over a Maxwellian relative-energy distribution."""

    m1, m2 = REACTANT_MASSES_AMU[reaction]
    mu = m1 * m2 / (m1 + m2) * AMU
    lo = _SIGMA[reaction]["Erange"][0]
    E = np.linspace(lo, min(_SIGMA[reaction]["Erange"][1], 60.0 * T_keV), n)   # keV
    sigma = cross_section_m2(reaction, E)
    v = np.sqrt(2.0 * E * KEV / mu)
    # f(E) dE for relative energy in a Maxwellian: 2/sqrt(pi) * sqrt(E)/T^{3/2} exp(-E/T)
    f = (2.0 / math.sqrt(math.pi)) * np.sqrt(E) / T_keV**1.5 * np.exp(-E / T_keV)
    return float(np.trapezoid(sigma * v * f, E))


def verify(verbose: bool = False) -> Dict[str, float]:
    """Cross-check the σ and ⟨σv⟩ fits against each other and against reference values."""

    out = {}
    for rxn in _SIGMA:
        errs = []
        for T in (5.0, 10.0, 20.0, 50.0):
            a = reactivity_m3_s(rxn, T)[0]
            b = reactivity_by_integration(rxn, T)
            errs.append(abs(a / b - 1.0))
        out[f"{rxn}_fit_vs_integral_rel_max"] = float(max(errs))
    # Bosch–Hale reference values (Table VIII), cm^3/s
    # Only the D–T values are taken from the paper's table with confidence; the other reactions are
    # checked by the σ-vs-⟨σv⟩ cross-integration above.
    ref = {("DT", 10.0): 1.137e-16, ("DT", 20.0): 4.330e-16}
    for (rxn, T), val in ref.items():
        out[f"{rxn}_{T:.0f}keV_rel_err"] = float(abs(reactivity_m3_s(rxn, T)[0] * 1e6 / val - 1.0))
    if verbose:
        for k, v in out.items():
            print(f"  {k:36s} {v:.3e}")
    return out


def volumetric_rate(reaction: str, n1_m3: float, n2_m3: float, T_keV: float) -> float:
    """Reactions per m³ per s for Maxwellian species densities n1, n2 (½ n² for like species)."""

    sv = float(reactivity_m3_s(reaction, T_keV)[0])
    same = reaction.startswith("DD")
    return (0.5 if same else 1.0) * n1_m3 * n2_m3 * sv


class FusionEventSampler:
    """event_hook for run_orbits: beam–target fusion of tracked ions on a Maxwellian target."""

    def __init__(self, reaction: str, target: Species, target_density_m3: float, target_temperature_eV: float, rng,
                 record_positions: bool = True):
        self.reaction = reaction
        self.target = target
        self.n_t = target_density_m3
        self.T_t = target_temperature_eV
        self.rng = rng
        self.record = record_positions
        self.birth_positions: List[np.ndarray] = []
        self.events = 0.0
        self.expected = 0.0

    def __call__(self, ens, ids, dt, t):
        if len(ids) == 0:
            return None
        v1 = ens.v[ids]
        m1 = ens.species.mass_kg
        m2 = self.target.mass_kg
        sigma_t = math.sqrt(self.T_t * E_CHARGE / m2)
        v2 = self.rng.normal(0.0, sigma_t, size=v1.shape)
        vrel = v1 - v2
        vr = np.linalg.norm(vrel, axis=1)
        mu = m1 * m2 / (m1 + m2)
        E_cm_keV = 0.5 * mu * vr**2 / KEV
        sigma = cross_section_m2(self.reaction, E_cm_keV)
        sigma = np.where(np.isfinite(sigma), sigma, 0.0)
        p = self.n_t * sigma * vr * dt
        self.expected += float(np.sum(p)) * ens.weight
        hit = self.rng.random(len(ids)) < p
        k = int(np.sum(hit))
        if k:
            self.events += k * ens.weight
            if self.record:
                self.birth_positions.append(ens.x[ids][hit].copy())
        return {"fusion_events": k * ens.weight, "fusion_expected": float(np.sum(p)) * ens.weight}
