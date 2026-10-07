"""Scale presets: benchtop, reactor, dimensionless.

Each preset fixes the vessel size, conductor, current, gas/plasma background,
the test-particle species and energies, and the integration window.  Field
grids are computed once per geometry at unit scale and rescaled exactly
(B ∝ I/L), so the three presets share one set of grids.

Benchtop  — the 2025 apartment prototype class: ~15 cm vessel radius, 1 mm
            copper conductor (auto-thinned when a winding's clearance demands
            it), 1 kA pulsed current, glow-discharge densities (10¹⁶ m⁻³,
            T_e ≈ 3 eV), D⁺ test ions at 10 eV, electrons at 3 eV.
Reactor   — the Fusion-lab normalisation: 1 m vessel, 5 T-class fields
            (superconducting-cable current, 1 cm conductor), n = 10²⁰ m⁻³,
            T = 15 keV, D⁺ at 15 keV, alphas at 3.5 MeV.
Dimensionless — the Codex trials' normalisation: unit ball, B_rms = 1 in the
            region of interest, speed 0.18, time in gyro-units; retention
            statistics directly comparable to the 162-run factorial.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Dict, Optional

from .constants import ALPHA, AMU, DEUTERON, ELECTRON, E_CHARGE, MU0, TRITON, Species

# Dimensionless units: choose a Species with m = 1, q = 1 in "code units" by
# building an ad-hoc species and interpreting SI arrays as code units.
CODE_UNIT_SPECIES = Species("code+", 1.0, 1.0)


@dataclass
class ScalePreset:
    name: str
    vessel_radius_m: float           # wall sphere radius
    unit_ball_m: float               # scale of the winding generation ball (wall = 0.82 × this)
    wire_radius_m: float
    current_A: float
    ramp_time_s: Optional[float]     # current ramp for induced-E runs (None = static)
    background: Dict[str, dict]      # name -> {species, density_m3, temperature_eV}
    test_species: Dict[str, dict]    # label -> {species, energy_eV, n}
    t_max_s: float
    dt_s: float
    fuel: Optional[dict] = None      # fusion target description
    notes: str = ""
    units: str = "SI"

    def describe(self) -> str:
        return (f"{self.name}: R_wall={self.vessel_radius_m:.3g} m, a_wire={self.wire_radius_m:.3g} m, I={self.current_A:.3g} A, "
                f"t_max={self.t_max_s:.3g} s, dt={self.dt_s:.3g} s")


def benchtop() -> ScalePreset:
    return ScalePreset(
        name="benchtop",
        vessel_radius_m=0.82 * 0.18, unit_ball_m=0.18, wire_radius_m=1.0e-3, current_A=1000.0, ramp_time_s=1.0e-3,
        background={
            "electrons": {"species": ELECTRON, "density_m3": 1e16, "temperature_eV": 3.0},
            "argon_ions": {"species": Species("Ar+", 39.948 * AMU, E_CHARGE), "density_m3": 1e16, "temperature_eV": 0.1},
        },
        test_species={
            "D+ 10 eV": {"species": DEUTERON, "energy_eV": 10.0, "n": 256},
            "e- 3 eV": {"species": ELECTRON, "energy_eV": 3.0, "n": 128},
        },
        t_max_s=2.0e-5, dt_s=1.0e-7,
        fuel={"reaction": "DD_nHe3", "target": DEUTERON, "density_m3": 1e16, "temperature_eV": 0.1},
        notes="Pulsed 1 kA in copper; a glow discharge at ~4 Pa. Fusion rates here are astronomically small and are reported only as a ledger sanity check.",
    )


def reactor() -> ScalePreset:
    # 1 m unit ball; current chosen so the reference solenoid gives ~5 T on axis with 12 turns over 1.1 m: B ≈ μ0 N I / L
    I = 5.0 * 1.1 / (MU0 * 12)
    return ScalePreset(
        name="reactor",
        vessel_radius_m=0.82, unit_ball_m=1.0, wire_radius_m=1.0e-2, current_A=I, ramp_time_s=None,
        background={
            "electrons": {"species": ELECTRON, "density_m3": 1e20, "temperature_eV": 15e3},
            "deuterons": {"species": DEUTERON, "density_m3": 5e19, "temperature_eV": 15e3},
            "tritons": {"species": TRITON, "density_m3": 5e19, "temperature_eV": 15e3},
        },
        test_species={
            "D+ 15 keV": {"species": DEUTERON, "energy_eV": 15e3, "n": 256},
            "alpha 3.5 MeV": {"species": ALPHA, "energy_eV": 3.5e6, "n": 128},
        },
        t_max_s=2.0e-5, dt_s=2.0e-9,
        fuel={"reaction": "DT", "target": TRITON, "density_m3": 5e19, "temperature_eV": 15e3},
        notes=f"I = {I:.3g} A per conductor so the 12-turn reference solenoid reaches ~5 T; Chad-core windings see the same current.",
    )


def dimensionless() -> ScalePreset:
    """Codex normalisation: unit ball, |B|_rms(ROI) = 1, speed 0.18, q = m = 1."""

    return ScalePreset(
        name="dimensionless",
        vessel_radius_m=0.82, unit_ball_m=1.0, wire_radius_m=0.005, current_A=1.0, ramp_time_s=None,
        background={},
        test_species={"code particles": {"species": CODE_UNIT_SPECIES, "energy_eV": None, "n": 128}},
        t_max_s=19.2, dt_s=0.008, fuel=None, units="code",
        notes="Speed 0.18, wall 0.82, t_max 19.2, dt 0.008 as in hopf-factorial-data.json; current normalised per run to B_rms = 8 in the ROI (their reference).",
    )


PRESETS = {"benchtop": benchtop, "reactor": reactor, "dimensionless": dimensionless}
