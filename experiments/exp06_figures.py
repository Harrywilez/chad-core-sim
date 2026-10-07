"""Figures for REPORT.md from the saved results."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ccsim.runner import unit_catalogue

R = ROOT / "results"
F = ROOT / "figures"
F.mkdir(exist_ok=True)
cat = unit_catalogue()
short = {w.name: w.name.replace(" (10 circuits)", "").replace("recursive-", "rec ").replace("codex ", "cx ").replace("hopf-drift ", "hd ") for w in cat}

# 1. geometry gallery
fig, axes = plt.subplots(3, 5, figsize=(16, 9.5), subplot_kw={"projection": "3d"})
for ax, w in zip(axes.ravel(), cat):
    p = w.active_points
    ax.plot(p[:, 0], p[:, 1], p[:, 2], lw=0.5, color="tab:blue")
    ax.set_title(short[w.name], fontsize=8)
    ax.set_xlim(-1, 1); ax.set_ylim(-1, 1); ax.set_zlim(-1, 1)
    ax.set_axis_off()
for ax in axes.ravel()[len(cat):]:
    ax.set_axis_off()
fig.suptitle("Chad-core configurations (active winding, unit ball; remote return not shown)", fontsize=11)
fig.tight_layout()
fig.savefig(F / "fig01_geometries.png", dpi=130)

# 2. topology + retention comparison (survival after 4 crossing times = equal path length; and end-of-window)
topo = json.load(open(R / "exp03_topology.json"))
T = json.load(open(R / "exp04_theory_tests.json"))
rows = {r["winding"]: r for r in T["T1"]["rows"]}
names = [w.name for w in cat]
pred = [topo[n]["predicted_adiabatic_retention"] for n in names]
fig, axes = plt.subplots(2, 1, figsize=(13, 9), sharex=True)
x = np.arange(len(names))
species = (("bench_e", "benchtop e⁻ 3 eV", "benchtop_e_adiabaticity", "tab:green"),
           ("reactor_D", "reactor D⁺ 15 keV", "reactor_D_adiabaticity", "tab:blue"),
           ("reactor_alpha", "reactor α 3.5 MeV", "reactor_alpha_adiabaticity", "tab:purple"),
           ("dimless", "dimensionless (Codex normalisation)", "dimensionless_adiabaticity", "tab:orange"),
           ("bench_D", "benchtop D⁺ 10 eV", "benchtop_D_adiabaticity", "tab:red"))
for ax, key, title in ((axes[0], "S4", "survival after 4 wall-transit times (equal path length for every species and scale)"),
                       (axes[1], "S8", "survival after 8 wall-transit times")):
    ax.bar(x - 0.3, pred, 0.18, label="adiabatic prediction from field-line mirror ratios", color="k")
    for i, (k, label, akey, color) in enumerate(species):
        vals = [rows[n][key].get(k) if rows[n][key].get(k) is not None else np.nan for n in names]
        adiab = np.nanmedian([rows[n][akey] for n in names if rows[n].get(akey) is not None])
        n_m = 128 if k in ("bench_e", "reactor_alpha", "dimless") else 256
        err = 1.96 * np.sqrt(np.array(vals) * (1 - np.array(vals)) / n_m)
        ax.errorbar(x - 0.15 + 0.1 * i, vals, yerr=err, fmt="o", ms=4, color=color, capsize=2, label=f"{label}  (ρ_L/L_B ≈ {adiab:.0e})")
    ax.set_ylabel("fraction still confined")
    ax.set_title(title, fontsize=10)
axes[0].legend(fontsize=7, ncol=2)
axes[1].set_xticks(x)
axes[1].set_xticklabels([short[n] for n in names], rotation=35, ha="right", fontsize=8)
fig.suptitle("T1 — field-line prediction vs full orbits: mirror trapping is real, then drifts empty the trap", fontsize=11)
fig.tight_layout()
fig.savefig(F / "fig02_prediction_vs_orbits.png", dpi=150)

# 3. adiabaticity sweep
th = R / "exp04_theory_tests.json"
if th.exists():
    T = json.load(open(th))
    sw = T["T2"]["sweep"]
    fig, ax = plt.subplots(figsize=(7, 4.2))
    a = [s["adiabaticity"] for s in sw]
    r = [s["retention"] for s in sw]
    lo = [s["wilson"][0] for s in sw]
    hi = [s["wilson"][1] for s in sw]
    ax.errorbar(a, r, yerr=[np.array(r) - np.array(lo), np.array(hi) - np.array(r)], fmt="o-", capsize=3, label="full orbits, equal path length")
    ax.axhline(T["T2"]["predicted_adiabatic"], color="k", ls="--", label="adiabatic prediction")
    codex_pt = next((s for s in sw if abs(s["speed"] - 0.18) < 1e-9), None)
    if codex_pt:
        ax.plot([codex_pt["adiabaticity"]], [codex_pt["retention"]], "r*", ms=14, label="Codex operating point (speed 0.18)")
    ax.set_xscale("log")
    ax.set_xlabel("adiabaticity ρ_L / L_B (median over ensemble)")
    ax.set_ylabel("retention")
    ax.set_title(f"T2 — {short[T['T2']['winding']]}: retention vs gyroradius", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(F / "fig03_adiabaticity_sweep.png", dpi=150)

    # 4. Codex comparison
    comp = T["T3"]
    if comp:
        fig, ax = plt.subplots(figsize=(7, 4))
        xx = np.arange(len(comp))
        ax.errorbar(xx - 0.12, [c["codex_retention"] for c in comp], yerr=[[c["codex_retention"] - c["codex_ci"][0] for c in comp], [c["codex_ci"][1] - c["codex_retention"] for c in comp]],
                    fmt="s", capsize=3, label="Codex (midpoint Biot–Savart, softened)")
        ax.errorbar(xx + 0.12, [c["ccsim_retention"] for c in comp], yerr=[[c["ccsim_retention"] - c["ccsim_wilson"][0] for c in comp], [c["ccsim_wilson"][1] - c["ccsim_retention"] for c in comp]],
                    fmt="o", capsize=3, label="ccsim (exact segments, finite conductor)")
        ax.set_xticks(xx)
        ax.set_xticklabels([c["variant"] for c in comp])
        ax.set_ylabel("retention (128 particles, t = 19.2)")
        ax.set_title("T3 — same normalisation, two independent implementations", fontsize=10)
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(F / "fig04_codex_comparison.png", dpi=150)

# 5. wall-flux maps
ch = R / "exp05_chemistry.json"
if ch.exists():
    C = json.load(open(ch))
    wf = C["wall_flux"]
    pick = [n for n in names if n in wf][:6]
    fig, axes = plt.subplots(2, 3, figsize=(13, 6))
    for ax, n in zip(axes.ravel(), pick):
        m = np.array(wf[n]["flux_map"])
        im = ax.imshow(m / max(m.mean(), 1e-300), origin="lower", aspect="auto", extent=[-180, 180, -1, 1], cmap="magma")
        ax.set_title(f"{short[n]}  peak/mean={wf[n]['peak_to_mean']:.1f}", fontsize=8)
        ax.set_xlabel("φ (deg)"); ax.set_ylabel("cos θ")
    fig.colorbar(im, ax=axes.ravel().tolist(), label="ion flux / mean")
    fig.suptitle("benchtop D⁺ wall-flux maps (ion bombardment pattern seen by a catalyst liner)", fontsize=10)
    fig.savefig(F / "fig05_wall_flux.png", dpi=140)

# 6. field magnitude slices for three windings
from ccsim.runner import build_grids
grids = build_grids(cat, ROOT / "cache", verbose=False)
fig, axes = plt.subplots(1, 4, figsize=(16, 4))
for ax, key in zip(axes, ["recursive-L3 (coil of coils, 12×6)", "codex Hopf-coordinate drift", "hopf-drift Shallow 180° (10 circuits)", "reference solenoid (12 turns)"]):
    g = grids[key]
    k = g.n // 2
    sl = g.Bmag[:, :, k]
    im = ax.imshow(np.log10(np.maximum(sl.T, 1e-9)), origin="lower", extent=[-g.half, g.half, -g.half, g.half], cmap="viridis", vmin=-7, vmax=-4.5)
    circ = plt.Circle((0, 0), 0.82, fill=False, color="w", lw=0.8, ls="--")
    ax.add_patch(circ)
    ax.set_title(short[key], fontsize=9)
    ax.set_xlabel("x"); ax.set_ylabel("y")
fig.colorbar(im, ax=axes.ravel().tolist(), label="log10 |B| (T per A, unit ball)")
fig.suptitle("|B| in the z = 0 plane (dashed: wall)", fontsize=10)
fig.savefig(F / "fig06_field_slices.png", dpi=130)
print("figures written")
