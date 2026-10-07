"""Orchestration: grids for the catalogue, particle screens per preset, ledgers."""

from __future__ import annotations

import json
import math
import os
import time
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np

from .collisions import Background, CoulombOperator
from .constants import E_CHARGE
from .fields import FieldGrid, biot_savart, return_sensitivity
from .fusion import FusionEventSampler, reactivity_m3_s
from .geometry import Winding, chadcore_catalogue, nearest_distance
from .particles import Ensemble, ScaledGridField, maxwellian_ensemble, monoenergetic_ensemble, run_orbits, sample_positions_in_vessel
from .presets import ScalePreset


def unit_catalogue(include_hopf_drift: bool = True) -> List[Winding]:
    """All configurations generated in the unit ball with 1 A and a 5 mm conductor."""

    return chadcore_catalogue(1.0, 0.005, 1.0, include_hopf_drift=include_hopf_drift)


def build_grids(catalogue: List[Winding], cache_dir: Path, n: int = 41, half_extent: float = 0.87,
                with_A: bool = True, verbose: bool = True) -> Dict[str, FieldGrid]:
    """Compute (or load) one unit-scale grid per configuration."""

    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    grids = {}
    for w in catalogue:
        key = "".join(c if c.isalnum() else "_" for c in w.name)[:60]
        path = cache_dir / f"grid_{key}_n{n}.npz"
        t0 = time.time()
        if path.exists():
            data = np.load(path)
            g = FieldGrid.__new__(FieldGrid)
            g.winding, g.half, g.n = w, half_extent, n
            g.axis = np.linspace(-half_extent, half_extent, n)
            g.dx = float(g.axis[1] - g.axis[0])
            xx, yy, zz = np.meshgrid(g.axis, g.axis, g.axis, indexing="ij")
            g.points = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
            g.B = data["B"]
            g.A = data["A"] if "A" in data.files else None
            g.Bmag = np.linalg.norm(g.B, axis=3)
            g.wire_distance = data["wire_distance"]
            src = "cache"
        else:
            g = FieldGrid(w, half_extent, n, with_A=with_A)
            np.savez_compressed(path, B=g.B, A=g.A if g.A is not None else np.zeros(0), wire_distance=g.wire_distance)
            src = "computed"
        grids[w.name] = g
        if verbose:
            print(f"  grid {w.name:55s} {src:8s} {time.time() - t0:5.1f}s", flush=True)
    return grids


def roi_mask(grid: FieldGrid, r_min: float = 0.20, r_max: float = 0.65, wire_clear: float = 0.08) -> np.ndarray:
    r = np.linalg.norm(grid.points, axis=1)
    return (r >= r_min) & (r <= r_max) & (grid.wire_distance.ravel() >= wire_clear)


def field_summary(w: Winding, grid: FieldGrid) -> dict:
    mask = roi_mask(grid)
    B = grid.B.reshape(-1, 3)[mask]
    brms = float(np.sqrt(np.mean(np.sum(B * B, axis=1))))
    rng = np.random.default_rng(0)
    sample = grid.points[mask][rng.choice(np.flatnonzero(mask).size, size=min(150, int(mask.sum())), replace=False)]
    # field-line "mirror ratio" proxy: B at wall shell / B in core, and gradient scale
    core = (np.linalg.norm(grid.points, axis=1) < 0.25) & (grid.wire_distance.ravel() > 0.05)
    shell = (np.abs(np.linalg.norm(grid.points, axis=1) - 0.75) < 0.05) & (grid.wire_distance.ravel() > 0.05)
    return {
        "brms_roi_per_amp_per_m": brms,
        "b_core_median": float(np.median(grid.Bmag.ravel()[core])) if core.any() else float("nan"),
        "b_shell_median": float(np.median(grid.Bmag.ravel()[shell])) if shell.any() else float("nan"),
        "return_sensitivity": return_sensitivity(w, sample),
        "interp_error": grid.interpolation_error(),
        "active_length": w.active_length_m,
        "total_length": w.length_m,
        "clearance": w.clearance_m(),
    }


def make_field(grid: FieldGrid, preset: ScalePreset, unit_winding: Winding, current_A: float, wire_radius_m: float):
    """Scaled field for a preset; the conductor radius is thinned to fit the clearance if needed."""

    scale = preset.unit_ball_m
    clearance = unit_winding.clearance_m() * scale
    a = min(wire_radius_m, 0.45 * clearance)
    f = ScaledGridField(grid, scale, current_A)
    f.wire_radius = a
    return f, a


def initial_positions(n: int, rng, unit_catalogue_points: np.ndarray, r_min=0.18, r_max=0.58, min_wire=0.035) -> np.ndarray:
    return sample_positions_in_vessel(n, r_min, r_max, rng, avoid_points=unit_catalogue_points, min_distance=min_wire)


def run_configuration(w: Winding, grid: FieldGrid, preset: ScalePreset, positions_unit: np.ndarray, rng,
                      species_label: Optional[str] = None, collisions: bool = False, fusion: bool = False,
                      current_override: Optional[float] = None, t_max_override: Optional[float] = None) -> dict:
    scale = preset.unit_ball_m
    current = preset.current_A if current_override is None else current_override
    field, a = make_field(grid, preset, w, current, preset.wire_radius_m)
    out = {"winding": w.name, "family": w.meta.get("family"), "preset": preset.name, "current_A": current, "wire_radius_m": a,
           "wire_thinned": a < preset.wire_radius_m, "runs": {}}
    scaled = w.scaled(scale).with_current(current).with_wire_radius(a)
    out["resistance_ohm"] = scaled.resistance_ohm()
    out["joule_power_W"] = scaled.joule_power_W()
    out["current_density_A_m2"] = current / (math.pi * a * a)
    out["B_max_T"] = field.max_B()
    mask = roi_mask(grid)
    out["B_rms_roi_T"] = float(abs(field.gain) * np.sqrt(np.mean(grid.Bmag.ravel()[mask] ** 2)))
    for label, spec in preset.test_species.items():
        if species_label and label != species_label:
            continue
        sp = spec["species"]
        n = spec["n"]
        pos = positions_unit[:n] * scale
        if preset.units == "code":
            d = rng.normal(size=(n, 3))
            d /= np.linalg.norm(d, axis=1)[:, None]
            ens = Ensemble(sp, pos, 0.18 * d)
            t_max, dt = preset.t_max_s, preset.dt_s
        else:
            ens = monoenergetic_ensemble(sp, n, spec["energy_eV"], pos, rng)
            v = math.sqrt(2 * spec["energy_eV"] * E_CHARGE / sp.mass_kg)
            crossings = 8.0 if sp.name.startswith("e") else 40.0
            t_max = t_max_override or crossings * preset.vessel_radius_m / v
            # coarse step from the typical field; sub-stepping resolves the strong field near conductors
            omega_ref = abs(sp.charge_C) / sp.mass_kg * max(out["B_rms_roi_T"], 1e-12)
            dt = min(preset.dt_s, 0.15 / omega_ref)
        coll = None
        if collisions and preset.background:
            coll = CoulombOperator([Background(b["species"], b["density_m3"], b["temperature_eV"]) for b in preset.background.values()], rng)
        hook = None
        sampler = None
        if fusion and preset.fuel and sp.name.startswith("D"):
            f = preset.fuel
            sampler = FusionEventSampler(f["reaction"], f["target"], f["density_m3"], f["temperature_eV"], rng)
            hook = sampler
        t0 = time.time()
        res = run_orbits(ens, field, t_max, dt, preset.vessel_radius_m, collision_operator=coll, event_hook=hook)
        rec = {
            "species": sp.name, "n": n, "t_max": t_max, "dt": dt, "retention": res.retention, "wilson": [res.wilson_lo, res.wilson_hi],
            "counts": res.counts, "median_loss_time": res.median_loss_time, "median_censored": res.median_censored,
            "mean_loss_time_restricted": res.mean_loss_time_restricted, "energy_drift_rel_max": res.energy_drift_rel_max,
            "mu_variation_median": res.mu_variation_median, "adiabaticity_median": res.adiabaticity_median,
            "substep_max": res.substep_max, "seconds": time.time() - t0, "collisions": bool(coll), "fusion": bool(hook),
            "loss_positions": res.loss_positions.tolist(), "loss_kinds": res.loss_kinds.tolist(), "loss_times": res.loss_times.tolist(),
        }
        if sampler is not None:
            confined_marker_seconds = float(np.sum(np.where(np.isfinite(ens.loss_time), ens.loss_time, t_max)))
            rec["fusion_events"] = sampler.events
            rec["fusion_expected"] = sampler.expected
            rec["fusion_rate_per_marker_s"] = sampler.expected / max(confined_marker_seconds, 1e-300)
            T_keV = preset.fuel["temperature_eV"] / 1e3
            rec["maxwellian_sv_m3_s"] = float(reactivity_m3_s(preset.fuel["reaction"], T_keV)[0])
            rec["maxwellian_rate_per_ion_s"] = rec["maxwellian_sv_m3_s"] * preset.fuel["density_m3"]
        out["runs"][label] = rec
    return out


def normalise_dimensionless_current(grid: FieldGrid, target_brms: float = 8.0) -> float:
    mask = roi_mask(grid)
    brms = float(np.sqrt(np.mean(grid.Bmag.ravel()[mask] ** 2)))
    return target_brms / brms
