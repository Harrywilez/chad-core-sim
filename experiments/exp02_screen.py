"""Particle motion in every Chad-core configuration at all three scales.

For each of the 13 windings (3 recursive, 3 Codex, 6 Hopf-drift, 1 reference
solenoid) and each preset:
  * field summary (B_rms in the ROI, core/shell field ratio, return-path
    sensitivity, grid interpolation error, conductor clearance, Joule power)
  * full-orbit ensembles of the preset's test species, same initial
    positions for every winding, finite-window retention with Wilson interval,
    loss classification, median loss time, energy/μ diagnostics
  * reactor scale: D⁺ runs repeated with Coulomb collisions and with the
    beam–target D–T fusion sampler; alpha (3.5 MeV) runs
  * dimensionless: normalised exactly as the Codex factorial (B_rms = 8 in the
    ROI, speed 0.18, dt 0.008, t_max 19.2) so retention can be compared with
    hopf-factorial-data.json from an independent implementation.
Usage: python3 exp02_screen.py [benchtop|reactor|dimensionless|all]
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ccsim.presets import PRESETS
from ccsim.runner import build_grids, field_summary, initial_positions, normalise_dimensionless_current, run_configuration, unit_catalogue

which = sys.argv[1] if len(sys.argv) > 1 else "all"
names = ["benchtop", "reactor", "dimensionless"] if which == "all" else [which]

cat = unit_catalogue()
grids = build_grids(cat, ROOT / "cache", n=41, with_A=True, verbose=False)
rng = np.random.default_rng(271828)
union = np.vstack([w.active_points[::4] for w in cat])
positions = initial_positions(256, rng, union)

results_dir = ROOT / "results"
results_dir.mkdir(exist_ok=True)
summary_path = results_dir / "exp02_field_summary.json"
if summary_path.exists():
    field_summaries = json.load(open(summary_path))
else:
    field_summaries = {}
    for w in cat:
        field_summaries[w.name] = field_summary(w, grids[w.name])
        print(f"field {w.name:55s} Brms/A/m={field_summaries[w.name]['brms_roi_per_amp_per_m']:.3e}  clearance={field_summaries[w.name]['clearance']:.4f}  "
              f"return-sens={field_summaries[w.name]['return_sensitivity']:.3f}  interp p95={field_summaries[w.name]['interp_error']['rel_p95']:.3e}", flush=True)
    json.dump(field_summaries, open(summary_path, "w"), indent=1)

for pname in names:
    preset = PRESETS[pname]()
    out_path = results_dir / f"exp02_{pname}.json"
    records = json.load(open(out_path)) if out_path.exists() else {}
    print(f"\n=== {preset.describe()}", flush=True)
    for w in cat:
        if w.name in records:
            continue
        g = grids[w.name]
        seed = np.random.default_rng(12345)
        t0 = time.time()
        if pname == "dimensionless":
            I = normalise_dimensionless_current(g, 8.0)
            rec = run_configuration(w, g, preset, positions, seed, current_override=I)
        elif pname == "reactor":
            rec = run_configuration(w, g, preset, positions, seed)
            rec["runs_collisional"] = run_configuration(w, g, preset, positions, np.random.default_rng(777), species_label="D+ 15 keV", collisions=True, fusion=True)["runs"]
        else:
            rec = run_configuration(w, g, preset, positions, seed)
        rec["field"] = field_summaries[w.name]
        records[w.name] = rec
        line = " | ".join(f"{k}: ret={v['retention']:.3f} [{v['wilson'][0]:.2f},{v['wilson'][1]:.2f}] wall={v['counts']['wall']} wire={v['counts']['wire']} sub={v['substep_max']}"
                         for k, v in rec["runs"].items())
        print(f"{w.name:55s} B_rms={rec['B_rms_roi_T']:.3g}  a={rec['wire_radius_m']:.2e}  P={rec['joule_power_W']:.3g} W  {line}  ({time.time()-t0:.0f}s)", flush=True)
        json.dump(records, open(out_path, "w"))
print("done")
