"""Chemistry coupled to the field: benchtop CO₂/H₂ discharge in a Chad-core winding.

1. Wall-flux map: where do the benchtop D⁺ (proxy for CO₂⁺/H₂⁺) markers land?
   The particle screen's loss positions are binned into flux per wall patch;
   the peak/mean ratio says how non-uniform the ion bombardment of a catalyst
   liner would be for each winding.
2. Gas-phase + electron-impact network (co2_h2_plasma_network) integrated at
   n_e = 1e16 m^-3, T_e = 3 eV, 4 Pa CO₂:H₂ = 1:4 for 10 ms — CO / O / H
   production; sensitivity to T_e.
3. Surface network on a catalyst patch fed by the wall flux (illustrative
   Sabatier barriers, exact Eckart tunnelling): CH₄ turnover vs patch flux,
   deterministic vs Gillespie on 500 sites.
4. Q-Surface import: the methane/Pt(111) network exported from qsurf, with
   every sech2 barrier converted to (V0, dE, ω‡) — verifies that the imported
   thermal rates reproduce qsurf's own numbers at 900 K.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ccsim.chemistry import WallFluxMap, co2_h2_plasma_network, load_qsurface_network, sabatier_surface_network
from ccsim.constants import E_CHARGE
from ccsim.presets import PRESETS

R = ROOT / "results"
out = {}
bench = json.load(open(R / "exp02_benchtop.json"))
p = PRESETS["benchtop"]()

# 1. wall flux maps: assume each marker represents a physical ion flux from a 1e16 m^-3 plasma of D+ at 10 eV
n_ion = 1e16
v_ion = math.sqrt(2 * 10.0 * E_CHARGE / (2.014 * 1.66053906660e-27))
vessel_volume = 4 / 3 * math.pi * p.vessel_radius_m**3
maps = {}
for name, rec in bench.items():
    run = rec["runs"]["D+ 10 eV"]
    n_markers = run["n"]
    weight = n_ion * vessel_volume / n_markers       # ions per marker
    fm = WallFluxMap.from_losses(np.array(run["loss_positions"]), np.array(run["loss_kinds"]), run["t_max"], weight, p.vessel_radius_m)
    maps[name] = {"total_ion_rate_s": fm.total_rate_s, "peak_to_mean": fm.peak_to_mean(),
                  "mean_flux_m2_s": float(fm.flux_m2_s.mean()), "mean_ML_per_s": float(fm.monolayers_per_s().mean()),
                  "peak_ML_per_s": float(fm.monolayers_per_s().max()), "flux_map": fm.flux_m2_s.tolist()}
    print(f"wall flux {name:55s} mean={maps[name]['mean_ML_per_s']:.3g} ML/s  peak/mean={maps[name]['peak_to_mean']:.2f}")
out["wall_flux"] = maps

# 2. gas-phase network
net = co2_h2_plasma_network()
gas = {}
for Te in (2.0, 3.0, 5.0):
    t, y, k = net.integrate({"e": 1e16, "CO2": 1e21, "H2": 4e21, "M": 5e21}, 1e-2, 600.0, Te, fixed=("e", "M"))
    gas[str(Te)] = {"final": {s: float(v[-1]) for s, v in y.items()}, "rate_coefficients": dict(zip([r.id for r in net.reactions], map(float, k))),
                    "CO2_conversion": float(1 - y["CO2"][-1] / 1e21)}
    print(f"gas Te={Te} eV: CO2 conversion after 10 ms = {gas[str(Te)]['CO2_conversion']:.3f}; CO={y['CO'][-1]:.3g} H={y['H'][-1]:.3g} H2O={y['H2O'][-1]:.3g}")
out["gas_phase"] = gas

# 3. surface network fed by the wall flux of the best and worst windings
surf = {}
ranked = sorted(maps.items(), key=lambda kv: kv[1]["mean_ML_per_s"])
for name, m in (ranked[0], ranked[-1]):
    for T in (500.0, 600.0, 700.0):
        sn = sabatier_surface_network(m["mean_ML_per_s"] * 1e-3, 0.3, 0.1, T)   # ×1e-3: only ~0.1% of arrivals are CO2/H2 ions in a D+ proxy run
        t, y, _ = sn.integrate({"*": 1.0, "CO2(g)": 1.0, "H2(g)": 1.0}, 1000.0, T, fixed=("CO2(g)", "H2(g)"))
        ch4_rate = float(np.gradient(y["CH4(g)"], t)[-1])
        surf[f"{name}|{T:.0f}K"] = {"CH4_turnover_per_site_s": ch4_rate, "coverages": {s: float(v[-1]) for s, v in y.items() if s.endswith("*")}}
        print(f"surface {name[:30]:30s} T={T:.0f} K: CH4 TOF={ch4_rate:.3e} /site/s  θ*={y['*'][-1]:.3f}")
rng = np.random.default_rng(0)
sn = sabatier_surface_network(ranked[-1][1]["mean_ML_per_s"] * 1e-3, 0.3, 0.1, 600.0)
tg, hist, ev = sn.gillespie({"*": 500, "CO2(g)": 1, "H2(g)": 1}, 1000.0, 600.0, rng, volume_m3=500, fixed=("CO2(g)", "H2(g)"))
t, y, _ = sn.integrate({"*": 1.0, "CO2(g)": 1.0, "H2(g)": 1.0}, 1000.0, 600.0, fixed=("CO2(g)", "H2(g)"))
surf["gillespie_vs_ode_CH4_final"] = {"gillespie_500_sites": float(hist["CH4(g)"][-1]), "ode_x500": float(y["CH4(g)"][-1] * 500), "events": ev}
print("Gillespie vs ODE CH4:", surf["gillespie_vs_ode_CH4_final"])
out["surface"] = surf

# 4. Q-Surface import
qpath = ROOT / "reference" / "methane_pt111.json"
if qpath.exists():
    qn = load_qsurface_network(str(qpath), flux_ml_per_s=0.05, sticking=0.0416)
    k = qn.rate_coefficients(900.0)
    imported = {r.id: float(kk) for r, kk in zip(qn.reactions, k)}
    qsurf_ref = {"ch3_dehydrogenation": 1.142e9, "ch2_dehydrogenation": 1.968e9, "ch_dehydrogenation": 1.7005e7, "h2_recombination": 2.2034e8}
    out["qsurface_import"] = {"imported_k_900K": imported, "qsurf_reference_k_900K": qsurf_ref,
                              "rel_diff": {kk: imported[kk] / v - 1 for kk, v in qsurf_ref.items() if kk in imported}}
    print("Q-Surface import rel diff vs qsurf at 900 K:", out["qsurface_import"]["rel_diff"])

json.dump(out, open(R / "exp05_chemistry.json", "w"), indent=1, default=float)
print("saved")
