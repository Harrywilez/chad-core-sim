"""Theory tests on the screen results.

T1  Adiabatic prediction vs full orbits.  exp03 predicts, from field-line
    mirror ratios alone, the isotropic trapped fraction of adiabatic particles.
    The reactor-scale D+ runs (ρ_L/L_B ~ 0.01) should match it winding by
    winding; the dimensionless runs (ρ_L/L_B ~ 0.1–0.3) should NOT — they sit
    in the non-adiabatic regime where "retention" is transient meandering.
T2  Retention vs adiabaticity.  For one winding at dimensionless scale, sweep
    the particle speed (ρ_L ∝ v) over two decades.  Prediction: retention → the
    adiabatic mirror value as ρ_L/L_B → 0 and → 0 (free streaming) as
    ρ_L/L_B → ∞, with a non-monotonic hump in between; the Codex operating
    point (speed 0.18) lies on that hump.
T3  Independent-code comparison.  Same six Hopf-drift variants, same
    normalisation (B_rms = 8, speed 0.18, dt 0.008, t_max 19.2, 10 circuits,
    r = 0.005): retention from exact-segment Biot–Savart (this code) vs the
    Codex midpoint-softened Biot–Savart.  Prediction: agreement within the
    Wilson intervals; residual differences trace to the wire model.
T4  Fusion ledger at reactor scale.  Beam–target D–T rate per confined marker
    vs the Maxwellian n⟨σv⟩ at 15 keV, and the total confined-ion fusion power
    the configuration would produce at n = 1e20 given its retention.
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

from ccsim.constants import E_CHARGE
from ccsim.fusion import REACTION_ENERGY_MEV, reactivity_m3_s
from ccsim.particles import Ensemble, run_orbits
from ccsim.presets import PRESETS
from ccsim.runner import build_grids, initial_positions, make_field, normalise_dimensionless_current, unit_catalogue

R = ROOT / "results"
topo = json.load(open(R / "exp03_topology.json"))
reactor = json.load(open(R / "exp02_reactor.json"))
dimless = json.load(open(R / "exp02_dimensionless.json"))
bench = json.load(open(R / "exp02_benchtop.json"))
out = {}

# ---------------------------------------------------------------- T1
def survival_at(rec, crossings, R_wall, v):
    """Fraction still confined after `crossings` wall-radius transit times."""
    if rec is None:
        return None
    t = crossings * R_wall / v
    lost = np.array(rec["loss_times"], dtype=float)
    n = rec["n"]
    return float(1.0 - np.sum(lost <= min(t, rec["t_max"])) / n)

import math as _m
from ccsim.constants import DEUTERON as _D, ALPHA as _A, ELECTRON as _E
v_D = _m.sqrt(2 * 15e3 * E_CHARGE / _D.mass_kg); v_a = _m.sqrt(2 * 3.5e6 * E_CHARGE / _A.mass_kg)
v_e = _m.sqrt(2 * 3.0 * E_CHARGE / _E.mass_kg); v_bD = _m.sqrt(2 * 10.0 * E_CHARGE / _D.mass_kg)
rows = []
for name, t in topo.items():
    rec_r = reactor.get(name, {}).get("runs", {}).get("D+ 15 keV")
    rec_a = reactor.get(name, {}).get("runs", {}).get("alpha 3.5 MeV")
    rec_d = dimless.get(name, {}).get("runs", {}).get("code particles")
    rec_b = bench.get(name, {}).get("runs", {}).get("D+ 10 eV")
    rows.append({
        "winding": name, "predicted_adiabatic": t["predicted_adiabatic_retention"], "frac_closed": t["frac_closed"],
        "mirror_ratio_median": t["mirror_ratio_median"],
        "reactor_D_retention": rec_r["retention"] if rec_r else None, "reactor_D_wilson": rec_r["wilson"] if rec_r else None,
        "reactor_D_adiabaticity": rec_r["adiabaticity_median"] if rec_r else None,
        "reactor_alpha_retention": rec_a["retention"] if rec_a else None, "reactor_alpha_adiabaticity": rec_a["adiabaticity_median"] if rec_a else None,
        "dimensionless_retention": rec_d["retention"] if rec_d else None, "dimensionless_adiabaticity": rec_d["adiabaticity_median"] if rec_d else None,
        "benchtop_D_retention": rec_b["retention"] if rec_b else None, "benchtop_D_adiabaticity": rec_b["adiabaticity_median"] if rec_b else None,
        "benchtop_e_retention": (bench.get(name, {}).get("runs", {}).get("e- 3 eV") or {}).get("retention"),
        "benchtop_e_adiabaticity": (bench.get(name, {}).get("runs", {}).get("e- 3 eV") or {}).get("adiabaticity_median"),
        # survival after 4 and 8 crossing times (equal path length across species/scales)
        "S4": {"reactor_D": survival_at(rec_r, 4, 0.82, v_D), "reactor_alpha": survival_at(rec_a, 4, 0.82, v_a),
               "bench_e": survival_at(bench.get(name, {}).get("runs", {}).get("e- 3 eV"), 4, 0.148, v_e),
               "bench_D": survival_at(rec_b, 4, 0.148, v_bD), "dimless": survival_at(rec_d, 4, 0.82, 0.18)},
        "S8": {"reactor_D": survival_at(rec_r, 8, 0.82, v_D), "reactor_alpha": survival_at(rec_a, 8, 0.82, v_a),
               "bench_e": survival_at(bench.get(name, {}).get("runs", {}).get("e- 3 eV"), 8, 0.148, v_e),
               "bench_D": survival_at(rec_b, 8, 0.148, v_bD), "dimless": survival_at(rec_d, 8, 0.82, 0.18)},
        "S20": {"reactor_D": survival_at(rec_r, 20, 0.82, v_D), "reactor_alpha": survival_at(rec_a, 20, 0.82, v_a)},
    })
pred = np.array([r["predicted_adiabatic"] for r in rows if r["reactor_D_retention"] is not None])
obs = np.array([r["reactor_D_retention"] for r in rows if r["reactor_D_retention"] is not None])
inside = [abs(r["reactor_D_retention"] - r["predicted_adiabatic"]) <= max(0.5 * (r["reactor_D_wilson"][1] - r["reactor_D_wilson"][0]), 0.03)
          for r in rows if r["reactor_D_retention"] is not None]
out["T1"] = {"rows": rows, "mean_abs_dev_reactor": float(np.mean(np.abs(pred - obs))) if len(pred) else None,
             "fraction_within_interval": float(np.mean(inside)) if inside else None,
             "corr_pred_vs_reactor": float(np.corrcoef(pred, obs)[0, 1]) if len(pred) > 2 and pred.std() > 0 and obs.std() > 0 else None}
print("T1 adiabatic prediction vs survival after 4 / 8 / 20 crossings:")
for r in rows:
    s4, s8, s20 = r["S4"], r["S8"], r["S20"]
    print(f"  {r['winding']:52s} pred={r['predicted_adiabatic']:.3f} | S4: e={s4['bench_e']:.2f} D={s4['reactor_D']:.2f} a={s4['reactor_alpha']:.2f} dl={s4['dimless']:.2f} bD={s4['bench_D']:.2f}"
          f" | S8: e={s8['bench_e']:.2f} D={s8['reactor_D']:.2f} a={s8['reactor_alpha']:.2f} | S20: D={s20['reactor_D']:.2f} a={s20['reactor_alpha']:.2f} | end: D={r['reactor_D_retention']:.3f}")
pred4 = np.array([r["predicted_adiabatic"] for r in rows]); s4e = np.array([r["S4"]["bench_e"] for r in rows]); s4D = np.array([r["S4"]["reactor_D"] for r in rows])
out["T1"]["S4_corr_pred_vs_bench_e"] = float(np.corrcoef(pred4, s4e)[0, 1]); out["T1"]["S4_corr_pred_vs_reactor_D"] = float(np.corrcoef(pred4, s4D)[0, 1])
out["T1"]["S4_mean_abs_dev_bench_e"] = float(np.mean(np.abs(pred4 - s4e))); out["T1"]["S4_mean_abs_dev_reactor_D"] = float(np.mean(np.abs(pred4 - s4D)))
print("  mean |pred-obs| =", out["T1"]["mean_abs_dev_reactor"], " within-interval fraction =", out["T1"]["fraction_within_interval"])

# ---------------------------------------------------------------- T2
cat = unit_catalogue()
grids = build_grids(cat, ROOT / "cache", verbose=False)
rng = np.random.default_rng(271828)
union = np.vstack([w.active_points[::4] for w in cat])
positions = initial_positions(256, rng, union)
target = next(w for w in cat if "Shallow 180" in w.name)
g = grids[target.name]
preset = PRESETS["dimensionless"]()
I = normalise_dimensionless_current(g, 8.0)
field, a = make_field(g, preset, target, I, preset.wire_radius_m)
sweep = []
for speed in (0.01, 0.02, 0.045, 0.09, 0.18, 0.36, 0.72, 1.5, 3.0):
    rr = np.random.default_rng(9)
    n = 128
    d = rr.normal(size=(n, 3))
    d /= np.linalg.norm(d, axis=1)[:, None]
    ens = Ensemble(preset.test_species["code particles"]["species"], positions[:n], speed * d)
    t_max = 19.2 * 0.18 / speed          # same path length for every speed
    dt = min(0.008, 0.15 / 8.0)
    t0 = time.time()
    res = run_orbits(ens, field, t_max, dt, 0.82)
    sweep.append({"speed": speed, "retention": res.retention, "wilson": [res.wilson_lo, res.wilson_hi], "adiabaticity": res.adiabaticity_median,
                  "mu_variation": res.mu_variation_median, "counts": res.counts, "seconds": time.time() - t0})
    print(f"T2 speed={speed:5.3f} rho/L={res.adiabaticity_median:.4f} retention={res.retention:.3f} [{res.wilson_lo:.2f},{res.wilson_hi:.2f}] mu-var={res.mu_variation_median:.2f} ({time.time()-t0:.0f}s)", flush=True)
out["T2"] = {"winding": target.name, "predicted_adiabatic": topo[target.name]["predicted_adiabatic_retention"], "sweep": sweep}

# ---------------------------------------------------------------- T3
codex_path = Path("/mnt/user-data/uploads/01a059d3-efc5-74e0-97fe-0d7ab433bd53/hopf-factorial-data.json")
if not codex_path.exists():
    codex_path = ROOT / "reference" / "hopf-factorial-data.json"
comp = []
if codex_path.exists():
    codex = json.load(open(codex_path))
    runs = {r["runId"]: r for r in codex["runs"]}
    for name, rec in dimless.items():
        if "hopf-drift" not in name:
            continue
        variant = rec.get("field", {}) and next((k for k in ("s64_150", "s64_180", "s64_210", "d48_150", "d48_180", "d48_210")
                                                if k in json.dumps(rec) or True), None)
        # map by label
        label = name.split("hopf-drift ")[1].split(" (")[0]
        key = {"Shallow 150°": "s64_150", "Shallow 180°": "s64_180", "Shallow 210°": "s64_210", "Deep 150°": "d48_150", "Deep 180°": "d48_180", "Deep 210°": "d48_210"}[label]
        cr = runs.get(f"{key}-n10-r0.005-p1.0")
        mine = rec["runs"]["code particles"]
        if cr:
            comp.append({"variant": key, "codex_retention": cr["retention"], "codex_ci": [cr["retentionLo"], cr["retentionHi"]],
                         "codex_fieldRms": cr["fieldRms"], "codex_median_loss": cr["medianLoss"],
                         "ccsim_retention": mine["retention"], "ccsim_wilson": mine["wilson"], "ccsim_median_loss": mine["median_loss_time"],
                         "overlap": not (mine["wilson"][1] < cr["retentionLo"] or mine["wilson"][0] > cr["retentionHi"])})
            print(f"T3 {key}: codex ret={cr['retention']:.3f} [{cr['retentionLo']:.2f},{cr['retentionHi']:.2f}] median={cr['medianLoss']:.2f} | "
                  f"ccsim ret={mine['retention']:.3f} [{mine['wilson'][0]:.2f},{mine['wilson'][1]:.2f}] median={mine['median_loss_time']:.2f}  overlap={comp[-1]['overlap']}")
out["T3"] = comp

# ---------------------------------------------------------------- T4
ledger = []
p = PRESETS["reactor"]()
sv = float(reactivity_m3_s("DT", 15.0)[0])
for name, rec in reactor.items():
    c = rec.get("runs_collisional", {}).get("D+ 15 keV")
    if not c:
        continue
    n_fuel = p.fuel["density_m3"]
    vol = 4 / 3 * math.pi * p.vessel_radius_m**3
    maxwellian_power = 0.5e20 * 0.5e20 * sv * vol * REACTION_ENERGY_MEV["DT"] * 1e6 * E_CHARGE    # n_D n_T <sv> V E_f
    ledger.append({"winding": name, "retention_collisional": c["retention"], "wilson": c["wilson"],
                   "fusion_rate_per_marker_s": c.get("fusion_rate_per_marker_s"), "maxwellian_rate_per_ion_s": c.get("maxwellian_rate_per_ion_s"),
                   "beam_target_over_maxwellian": (c.get("fusion_rate_per_marker_s") or 0) / max(c.get("maxwellian_rate_per_ion_s") or 1e-300, 1e-300),
                   "mean_confinement_time_s": c["mean_loss_time_restricted"], "energy_confinement_time_needed_s_for_Q1": None,
                   "maxwellian_volume_power_W_if_confined": maxwellian_power})
    print(f"T4 {name:55s} ret(coll)={c['retention']:.3f} beam-target/Maxwellian={ledger[-1]['beam_target_over_maxwellian']:.2f} tau_mean={c['mean_loss_time_restricted']:.2e} s")
out["T4"] = {"ledger": ledger, "sv_DT_15keV_m3_s": sv, "note": "Lawson-type: at n=1e20, T=15 keV the D-T triple product needs tau_E ~ 3 s; the confinement times here are microseconds."}

json.dump(out, open(R / "exp04_theory_tests.json", "w"), indent=1, default=float)
print("saved")
