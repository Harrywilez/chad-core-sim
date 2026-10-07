"""Rung 3 (estimate): collective loss of a plasma in a toroidal field without transform.

In a pure toroidal field the ∇B/curvature drift v_d = (v∥² + v⊥²/2)/(Ω R) is vertical and
charge-dependent: ions drift one way, electrons the other.  The plasma polarizes on the
ion plasma period, and the resulting vertical E drives an E×B drift that moves the whole
column outward at about the ion drift speed (the charge separation cannot be shorted
along the field lines because they never connect top to bottom).  Loss time
    τ_pol ≈ a_plasma / v_d,     v_d ≈ 2 T / (e B R)   (thermal, T in eV)
A rotational transform ι shorts the separation along the field line if the connection
length L_c = 2πR/ι is short compared to the distance the column would drift in a
parallel transit: the shorting works when  ι ≳ v_d / v_th,i · (2πR / a)  — i.e. even a
small ι (0.05–0.3) is enough for thermal ions, after which the residual loss is
neoclassical diffusion rather than bulk drift.  This module just does that bookkeeping.
"""

from __future__ import annotations

import math

from .constants import E_CHARGE


def polarization_loss(T_eV: float, B_T: float, R_m: float, a_m: float, mass_kg: float, iota: float = 0.0,
                      shorting_species_mass_kg: float = 9.1093837e-31, iota_floor: float = 0.02) -> dict:
    """The charge separation is shorted along the field by the *electrons* (parallel current),
    so the shorting criterion compares the ion drift with the electron thermal speed over half
    a poloidal turn, L_c = πR/ι:  ι ≳ (v_d,i / v_th,e) · πR/a — plus a floor of iota_floor, below
    which the transform is in the noise of the field-line analysis."""
    v_th = math.sqrt(2 * T_eV * E_CHARGE / mass_kg)
    v_th_short = math.sqrt(2 * T_eV * E_CHARGE / shorting_species_mass_kg)
    v_d = 2 * T_eV / (B_T * R_m)                       # m/s  (2T/(eBR): thermal ∇B + curvature drift)
    tau_pol = a_m / v_d
    transit = R_m / v_th
    iota_needed = max(iota_floor, (v_d / v_th_short) * (math.pi * R_m / a_m))
    L_c = 2 * math.pi * R_m / iota if iota > 0 else float("inf")
    return {"v_drift_m_s": v_d, "tau_pol_s": tau_pol, "tau_pol_transits": tau_pol / transit,
            "iota_needed_to_short": iota_needed, "iota": iota, "connection_length_m": L_c,
            "shorted": bool(iota >= iota_needed), "v_th_m_s": v_th, "v_th_shorting_m_s": v_th_short}
