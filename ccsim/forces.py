"""Rung 4: electromagnetic forces on the conductors and the engineering ledger.

For every segment of every circuit: B from *all* segments (own straight segment
contributes nothing at its own midpoint; neighbouring segments use the uniform-current
interior model), then f = I t̂ × B  [N/m].  Reported per circuit and per crossing.

Stress proxies (order-of-magnitude, for a round conductor of radius a):
  * transverse load per length f  → bending/contact stress ~ f / (2a)   [Pa]
  * magnetic pressure at the conductor surface  B_s² / (2 μ0)          [Pa]
  * current density  J = I / (π a²)  against conductor limits
"""

from __future__ import annotations

import math

import numpy as np

from .constants import MU0
from .fastfield import SegmentSource

CONDUCTOR_LIMITS_A_M2 = {
    "copper, steady (water-cooled)": 2.0e7,
    "copper, pulsed (ms)": 1.0e9,
    "NbTi at 5 T, 4.2 K (engineering)": 3.0e8,
    "Nb3Sn at 12 T, 4.2 K (engineering)": 5.0e8,
    "REBCO tape at 20 T, 20 K (engineering)": 1.0e9,
}
YIELD_STRESS_PA = {"annealed copper": 7.0e7, "hard copper": 3.0e8, "316 stainless": 2.5e8, "Inconel 718": 1.0e9}


def segment_forces(winding, current_scale: float = 1.0) -> dict:
    src = SegmentSource(winding, current_scale)
    mid = src.p1 + 0.5 * src.dl
    B = src.B(mid)
    f = np.cross(src.that, B) * src.I[:, None]          # N/m
    fmag = np.linalg.norm(f, axis=1)
    Bs = np.linalg.norm(B, axis=1)
    # per-circuit bookkeeping
    parts = winding.windings if hasattr(winding, "windings") else [winding]
    counts = [len(w.points) - 1 for w in parts]
    idx = np.repeat(np.arange(len(parts)), counts)[: len(fmag)]
    per_circuit = []
    for k, w in enumerate(parts):
        m = idx == k
        net = np.sum(f[m] * src.L[m, None], axis=0)
        per_circuit.append({"name": w.name, "current_A": float(w.current_A * current_scale), "f_max_N_m": float(fmag[m].max()),
                            "f_mean_N_m": float(fmag[m].mean()), "B_surface_max_T": float(Bs[m].max()), "net_force_N": net.tolist()})
    a = float(np.min(src.a))
    I = float(np.max(np.abs(src.I)))
    return {
        "a_m": a, "I_A": I, "J_A_m2": I / (math.pi * a * a),
        "f_max_N_m": float(fmag.max()), "f_p99_N_m": float(np.percentile(fmag, 99)), "f_median_N_m": float(np.median(fmag)),
        "B_surface_max_T": float(Bs.max()), "magnetic_pressure_max_Pa": float(Bs.max() ** 2 / (2 * MU0)),
        "contact_stress_max_Pa": float(fmag.max() / (2 * a)),
        "per_circuit": per_circuit,
        "conductor_limits": {k: {"J_limit_A_m2": v, "ratio": I / (math.pi * a * a) / v} for k, v in CONDUCTOR_LIMITS_A_M2.items()},
        "yield_ratios": {k: float(fmag.max() / (2 * a) / v) for k, v in YIELD_STRESS_PA.items()},
        "segments": {"mid": mid, "f": f, "Bs": Bs},
    }


def crossing_force_estimate(I: float, d: float) -> float:
    """Force per unit length between two long straight antiparallel/perpendicular strands at
    separation d — the local scale of the crossing loads:  μ0 I² / (2π d)."""
    return MU0 * I * I / (2 * math.pi * d)
