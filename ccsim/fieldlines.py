"""Field-line topology and adiabatic confinement prediction.

Trace field lines from seeds in the region of interest in both directions
(RK4 on the unit tangent B/|B|), classify where they end (wall, conductor,
closed/returning, max length) and record the field strength along them.

For an adiabatic particle (ρ_L ≪ L_B) the motion is a bounce along the field
line between mirror points, so the confinement of an isotropic ensemble born
on a line is predicted by the mirror ratio in each direction:
    trapped  ⇔  sin²θ > B_seed / B_m,   B_m = min(max_+ B, max_− B)  (before the wall)
    P_trapped(isotropic) = sqrt(1 − B_seed / B_m)         (0 if a direction has no mirror)
This is the null hypothesis against which the full-orbit runs are tested:
a winding "confines" in the adiabatic limit only through mirror trapping on
open lines or through closed flux surfaces — no vacuum-field topology is
exempt from that.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Dict, List

import numpy as np

from .fields import FieldGrid


@dataclass
class LineResult:
    seed: np.ndarray
    end_plus: str
    end_minus: str
    length_plus: float
    length_minus: float
    bmax_plus: float
    bmax_minus: float
    b_seed: float
    closed: bool
    b_along_plus: np.ndarray
    s_along_plus: np.ndarray

    @property
    def mirror_ratio(self) -> float:
        if self.end_plus == "wall" and self.end_minus == "wall":
            return min(self.bmax_plus, self.bmax_minus) / self.b_seed
        if self.closed:
            return math.inf
        return min(self.bmax_plus, self.bmax_minus) / self.b_seed

    @property
    def predicted_trapped_fraction(self) -> float:
        if self.closed:
            return 1.0
        R = self.mirror_ratio
        return math.sqrt(max(0.0, 1.0 - 1.0 / R)) if R > 1.0 else 0.0


def trace_line(grid: FieldGrid, seed: np.ndarray, direction: float, wall_radius: float, step: float,
               max_length: float, wire_clear: float):
    p = np.asarray(seed, dtype=float).copy()
    s = 0.0
    bs = [float(grid.Bmag_at(p[None, :])[0])]
    ss = [0.0]
    end = "length"
    while s < max_length:
        def tangent(q):
            b = grid.B_at(q[None, :])[0]
            n = np.linalg.norm(b)
            return direction * b / n if n > 0 else np.zeros(3)
        k1 = tangent(p)
        k2 = tangent(p + 0.5 * step * k1)
        k3 = tangent(p + 0.5 * step * k2)
        k4 = tangent(p + step * k3)
        p = p + (step / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        s += step
        if not np.all(np.isfinite(p)) or np.linalg.norm(k1) == 0:
            end = "numerical"
            break
        r = np.linalg.norm(p)
        if r >= wall_radius:
            end = "wall"
            break
        if not np.all(np.abs(p) < grid.half):
            end = "grid-exit"
            break
        if grid.wire_distance_at(p[None, :])[0] <= wire_clear:
            end = "wire"
            break
        bs.append(float(grid.Bmag_at(p[None, :])[0]))
        ss.append(s)
        if s > 4 * step and np.linalg.norm(p - seed) < 0.5 * step and s > 1.0 * step * 10:
            end = "closed"
            break
    return end, s, np.array(bs), np.array(ss), p


def analyse(grid: FieldGrid, seeds: np.ndarray, wall_radius: float = 0.82, step: float = 0.008,
            max_length: float = 60.0, wire_clear: float = 0.01) -> List[LineResult]:
    out = []
    for seed in seeds:
        e_p, L_p, b_p, s_p, _ = trace_line(grid, seed, +1.0, wall_radius, step, max_length, wire_clear)
        e_m, L_m, b_m, s_m, _ = trace_line(grid, seed, -1.0, wall_radius, step, max_length, wire_clear)
        closed = e_p == "closed" or e_m == "closed" or (e_p == "length" and e_m == "length")
        out.append(LineResult(seed, e_p, e_m, L_p, L_m, float(b_p.max()), float(b_m.max()), float(b_p[0]), closed, b_p, s_p))
    return out


def summarise(lines: List[LineResult]) -> Dict[str, float]:
    ends = [l.end_plus for l in lines] + [l.end_minus for l in lines]
    n = len(lines)
    ratios = np.array([l.mirror_ratio for l in lines])
    finite = ratios[np.isfinite(ratios)]
    pred = np.array([l.predicted_trapped_fraction for l in lines])
    return {
        "n_seeds": n,
        "frac_wall_both": float(np.mean([l.end_plus == "wall" and l.end_minus == "wall" for l in lines])),
        "frac_wire_any": float(np.mean([l.end_plus == "wire" or l.end_minus == "wire" for l in lines])),
        "frac_closed": float(np.mean([l.closed for l in lines])),
        "frac_ends_wall": float(np.mean([e == "wall" for e in ends])),
        "connection_length_median": float(np.median([l.length_plus + l.length_minus for l in lines])),
        "mirror_ratio_median": float(np.median(finite)) if finite.size else float("nan"),
        "mirror_ratio_p90": float(np.percentile(finite, 90)) if finite.size else float("nan"),
        "predicted_adiabatic_retention": float(np.mean(pred)),
    }
