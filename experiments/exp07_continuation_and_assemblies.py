"""Hopf continuation ('keep revolving until it is a donut') and multi-core assemblies.

New configurations (unit ball):
  * hopf continued: sweep 720° / 48 circuits, 1440° / 96 circuits (η 0.82→0.64), 1440° / 96 (η→0.50)
  * hopf torus: η = 0.70 with 24 circuits × 1 revolution, 36 × 2
  * assemblies of the precessing-loop core in every layout (pairs, triad, tetra, octa cusp/aligned, cube)
  * assemblies of the Hopf-drift Shallow 180° core and of the reference solenoid in the key layouts
For each: field grid (33³), field-line topology + adiabatic prediction, dimensionless retention
(Codex normalisation) and reactor-scale D⁺ retention with survival at 4/8 transits.
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ccsim.constants import DEUTERON, E_CHARGE
from ccsim.fieldlines import analyse, summarise
from ccsim.fields import FieldGrid
from ccsim.geometry import (ASSEMBLY_LAYOUTS, MultiWinding, Winding, assembly_windings, close_with_return, codex_toroidal,
                            hopf_continued, hopf_drift, hopf_torus, solenoid)
from ccsim.presets import PRESETS
from ccsim.runner import initial_positions, normalise_dimensionless_current, run_configuration

R = ROOT / "results"
C = ROOT / "cache"
which = sys.argv[1] if len(sys.argv) > 1 else "all"


def single(name, pts, meta):
    closed, n_active = close_with_return(pts, 3.8)
    return Winding(name, closed, 1.0, 0.005, closed=True, active_count=n_active, meta=meta)


configs = []
if which in ("all", "hopf"):
    configs += [
        single("hopf continued 720° / 48 circuits", hopf_continued(0.82, 0.64, 720, 48), {"family": "hopf-continued", "sweep": 720, "circuits": 48, "eta1": 0.64}),
        single("hopf continued 1440° / 96 circuits", hopf_continued(0.82, 0.64, 1440, 96), {"family": "hopf-continued", "sweep": 1440, "circuits": 96, "eta1": 0.64}),
        single("hopf continued 1440° / 96, η→0.50", hopf_continued(0.82, 0.50, 1440, 96), {"family": "hopf-continued", "sweep": 1440, "circuits": 96, "eta1": 0.50}),
        single("hopf torus η=0.70, 24×1", hopf_torus(0.70, 24, 1.0), {"family": "hopf-torus", "eta": 0.70, "circuits": 24, "revolutions": 1}),
        single("hopf torus η=0.70, 36×2", hopf_torus(0.70, 36, 2.0), {"family": "hopf-torus", "eta": 0.70, "circuits": 36, "revolutions": 2}),
    ]
if which in ("all", "assemblies", "assemblies-precess", "assemblies-other"):
    cores = {
        "precess": codex_toroidal("precess", 0.22, 24.0, 3.0),
        "hopf-drift s64_180": hopf_drift("s64_180", 10),
        "solenoid": solenoid(0.45, 1.1, 12, n_per_turn=60),
    }
    key_layouts = ["pair-axial-same", "pair-axial-opposed", "pair-facing-flipped", "pair-side-crossed", "triad-120", "tetra-4", "octa-6-cusp", "cube-8-checker"]
    for cname, core in cores.items():
        if which == "assemblies-precess" and cname != "precess":
            continue
        if which == "assemblies-other" and cname == "precess":
            continue
        layouts = list(ASSEMBLY_LAYOUTS) if cname == "precess" else key_layouts
        for lay in layouts:
            if lay == "single":
                continue
            ws = assembly_windings(core, lay, 1.0, 0.005, 1.0, cname)
            configs.append(MultiWinding(f"{cname} × {lay}", ws, {"family": "assembly", "core": cname, "layout": lay, "cores": len(ws)}))

out_path = R / (f"exp07_{which}.json" if which != "all" else "exp07_continuation_and_assemblies.json")
results = json.load(open(out_path)) if out_path.exists() else {}
rng = np.random.default_rng(271828)
positions = initial_positions(256, rng, np.zeros((1, 3)), min_wire=0.0)
v_D = math.sqrt(2 * 15e3 * E_CHARGE / DEUTERON.mass_kg)

for w in configs:
    if w.name in results:
        continue
    t0 = time.time()
    key = "".join(c if c.isalnum() else "_" for c in w.name)[:60]
    path = C / f"grid_{key}_n33.npz"
    if path.exists():
        d = np.load(path)
        g = FieldGrid.__new__(FieldGrid)
        g.winding, g.half, g.n = w, 0.87, 33
        g.axis = np.linspace(-0.87, 0.87, 33); g.dx = float(g.axis[1] - g.axis[0])
        xx, yy, zz = np.meshgrid(g.axis, g.axis, g.axis, indexing="ij")
        g.points = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
        g.B = d["B"]; g.A = None; g.Bmag = np.linalg.norm(g.B, axis=3); g.wire_distance = d["wire_distance"]
    else:
        g = FieldGrid(w, 0.87, 33, with_A=False)
        np.savez_compressed(path, B=g.B, wire_distance=g.wire_distance)
    t_grid = time.time() - t0
    # field-line topology
    seeds = positions[:96]
    ok = g.wire_distance_at(seeds) > 0.03
    lines = analyse(g, seeds[ok], wall_radius=0.82, step=0.006, max_length=40.0, wire_clear=0.008)
    topo = summarise(lines)
    # dimensionless run
    pd = PRESETS["dimensionless"]()
    I = normalise_dimensionless_current(g, 8.0)
    rec_d = run_configuration(w, g, pd, positions, np.random.default_rng(12345), current_override=I)
    # reactor D+
    pr = PRESETS["reactor"]()
    rec_r = run_configuration(w, g, pr, positions, np.random.default_rng(12345), species_label="D+ 15 keV")
    rr = rec_r["runs"]["D+ 15 keV"]
    lost = np.array(rr["loss_times"], dtype=float)
    S = {f"S{c}": float(1 - np.sum(lost <= c * 0.82 / v_D) / rr["n"]) for c in (4, 8, 20)}
    for rec in (rec_d, rec_r):
        for run in rec["runs"].values():
            run.pop("loss_positions", None); run.pop("loss_times", None); run.pop("loss_kinds", None)
    results[w.name] = {"meta": w.meta, "clearance": w.clearance_m(), "length": w.length_m, "active_length": w.active_length_m,
                       "topology": {k: v for k, v in topo.items() if not k.startswith("per_seed")},
                       "dimensionless": rec_d, "reactor": rec_r, "reactor_survival": S, "grid_seconds": t_grid, "seconds": time.time() - t0}
    print(f"{w.name:45s} clr={w.clearance_m():.4f} wall={topo['frac_wall_both']:.2f} closed={topo['frac_closed']:.2f} Rmir p90={topo['mirror_ratio_p90']:.2f} "
          f"pred={topo['predicted_adiabatic_retention']:.3f} | dimless ret={rec_d['runs']['code particles']['retention']:.3f} | reactor D+ S4={S['S4']:.2f} S8={S['S8']:.2f} S20={S['S20']:.2f} "
          f"B_rms={rec_r['B_rms_roi_T']:.2f} T ({time.time()-t0:.0f}s)", flush=True)
    json.dump(results, open(out_path, "w"), default=float)
print("done")
