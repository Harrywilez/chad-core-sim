"""Physical constants (CODATA 2018 / SI 2019 exact values) and particle species."""

from __future__ import annotations

from dataclasses import dataclass

MU0 = 1.25663706212e-6          # N A^-2  (2018 CODATA)
EPS0 = 8.8541878128e-12         # F m^-1
C_LIGHT = 299792458.0           # m s^-1 (exact)
E_CHARGE = 1.602176634e-19      # C (exact)
K_B = 1.380649e-23              # J K^-1 (exact)
AMU = 1.66053906660e-27         # kg
M_E = 9.1093837015e-31          # kg
M_P = 1.67262192369e-27         # kg
EV = E_CHARGE                   # J per eV
KEV = 1.0e3 * EV
H_PLANCK = 6.62607015e-34       # J s (exact)
HBAR = H_PLANCK / (2.0 * 3.141592653589793)
N_A = 6.02214076e23

# Atomic masses in amu (nuclear masses would drop electron mass; the
# difference is < 0.1 % and irrelevant at the fidelity of a test-particle code).
MASS_AMU = {
    "e": M_E / AMU,
    "p": 1.007276,
    "D": 2.013553,
    "T": 3.015501,
    "He3": 3.014932,
    "He4": 4.001506,
    "C": 12.0,
    "O": 15.994915,
    "Ar": 39.948,
    "CO2": 43.98983,
    "CO": 27.994915,
    "H2": 2.01565,
    "H": 1.007825,
    "CH4": 16.0313,
    "H2O": 18.010565,
    "O2": 31.98983,
}


@dataclass(frozen=True)
class Species:
    """A charged (or neutral) particle species in SI units."""

    name: str
    mass_kg: float
    charge_C: float

    @property
    def mass_amu(self) -> float:
        return self.mass_kg / AMU

    @property
    def charge_number(self) -> float:
        return self.charge_C / E_CHARGE


def species(name: str, charge_number: int) -> Species:
    """Build a Species from the atomic-mass table, e.g. species('D', 1)."""

    if name not in MASS_AMU:
        raise KeyError(f"unknown species {name!r}; known: {sorted(MASS_AMU)}")
    return Species(f"{name}{'+' * charge_number if charge_number > 0 else '-' * (-charge_number) if charge_number < 0 else ''}",
                   MASS_AMU[name] * AMU, charge_number * E_CHARGE)


ELECTRON = Species("e-", M_E, -E_CHARGE)
DEUTERON = species("D", 1)
TRITON = species("T", 1)
HELION = species("He3", 2)
ALPHA = species("He4", 2)
PROTON = species("p", 1)
CO2_ION = species("CO2", 1)
ARGON_ION = species("Ar", 1)
