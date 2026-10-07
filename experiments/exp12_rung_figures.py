"""Figures for the rungs: fig12 Poincaré sections, fig13 transform vs modulation, fig14 GC survival, fig15 strand forces."""

from __future__ import annotations

import glob
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
R, F = ROOT / "results", ROOT / "figures"

surf = {}
for p in glob.glob(str(R / "exp10_*.json")):
    surf.update(json.load(open(p)))

# ---- fig12: Poincaré sections ----
picks = [("baseline (round torus)", "pure mirror pair, η 0.70 (toroidal field)"), ("helical-current ε=0.30, n=5", "chiral ε 0.30, n 5, η 0.70"),
         ("best: ε=0.70, n=4, η=0.60", "chiral ε 0.70, n 4, η 0.60"), ("best: ε=0.55, n=3, η=0.50", "chiral ε 0.55, n 3, η 0.50")]
plt.style.use("dark_background")
fig, axes = plt.subplots(1, 4, figsize=(16, 4.3))
col = {"surface": "#45b8a6", "island/chaotic": "#8a7bd6", "open": "#e05d57", "short": "#888"}
for ax, (key, title) in zip(axes, picks):
    if key not in surf:
        ax.set_title(f"{title}\n(missing)"); continue
    r = surf[key]
    for line, P in zip(r["per_line"], r["punctures"]):
        P = np.array(P)
        if len(P):
            ax.plot(P[:, 0], P[:, 1], ".", ms=2.2, color=col.get(line["kind"], "#888"), alpha=0.9)
    ax.plot(r["axis_R"], r["axis_z"], "+", color="w", ms=9, mew=1.2)
    ax.set_aspect("equal"); ax.set_facecolor("#0b1015"); ax.set_title(title, fontsize=9)
    ax.set_xlabel("R"); ax.set_ylabel("z")
    txt = f"surfaces {100*r['frac_surface']:.0f}% · ι {r['iota_axis']:+.2f}→{r['iota_edge']:+.2f}" if r["frac_surface"] > 0 else "no closed surfaces"
    ax.text(0.02, 0.02, txt, transform=ax.transAxes, fontsize=8, color="#ddd")
fig.suptitle("Poincaré sections at φ = 0 (exact Biot–Savart, 30–40 toroidal turns): teal = nested surface, purple = island/chaotic, red = open", fontsize=10)
fig.tight_layout(); fig.savefig(F / "fig12_poincare.png", dpi=140, facecolor="#0b1015"); plt.close(fig)
plt.style.use("default")

# ---- fig13: transform vs modulation ----
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
rows = []
for k, r in surf.items():
    m = r.get("meta", {})
    if m.get("family") == "hopf-mirror-chiral" and m.get("helical_mode") == "toroidal":
        rows.append((m["eps"], m["periods"], m["eta"], abs(r["iota_edge"]) if r["frac_surface"] > 0 else np.nan, r["outermost_surface_r"], r["frac_surface"]))
rows = np.array(rows, dtype=float)
for eta, mk in ((0.5, "s"), (0.6, "o"), (0.7, "^")):
    sel = rows[rows[:, 2] == eta]
    if len(sel):
        for n, c in ((3, "#d4884a"), (4, "#45b8a6"), (5, "#8a7bd6"), (6, "#e6b34f"), (7, "#7ad9ca")):
            s2 = sel[sel[:, 1] == n]
            if len(s2):
                o = np.argsort(s2[:, 0]); ax[0].plot(s2[o, 0], s2[o, 3], marker=mk, color=c, ls="-", label=f"η {eta}, n {n}")
ax[0].set_xlabel("helical density modulation ε"); ax[0].set_ylabel("|ι| at the last closed surface"); ax[0].grid(alpha=0.3); ax[0].legend(fontsize=7, ncol=2)
ax[0].set_title("rotational transform of the chiral meshed torus", fontsize=9)
ax[1].scatter(rows[:, 3], rows[:, 4], c=rows[:, 2], cmap="viridis", s=40)
for row in rows:
    ax[1].annotate(f"ε{row[0]:.2f} n{int(row[1])}", (row[3], row[4]), fontsize=6, alpha=0.7)
ax[1].set_xlabel("|ι| edge"); ax[1].set_ylabel("last closed surface radius (ball units)"); ax[1].grid(alpha=0.3)
ax[1].set_title("transform vs confined radius (colour = η: 0.5 dark → 0.7 bright)", fontsize=9)
fig.tight_layout(); fig.savefig(F / "fig13_transform_scan.png", dpi=140); plt.close(fig)

# ---- fig14: GC survival curves ----
tr = {}
for p in glob.glob(str(R / "exp11_*.json")):
    tr.update(json.load(open(p)))
if tr:
    fig, ax = plt.subplots(figsize=(7.5, 4.2))
    for name, r in tr.items():
        ax.plot(r["curve_t_transits"], r["curve_alive"], label=name)
    ax.set_xlabel("wall transits (15 keV D⁺, 1 m ball)"); ax.set_ylabel("fraction retained (guiding centre, exact field)"); ax.set_ylim(0, 1.02); ax.grid(alpha=0.3); ax.legend(fontsize=7)
    ax.set_title("rung 2: particles seeded inside the surface region — untwisted drifts out, chiral holds", fontsize=9)
    fig.tight_layout(); fig.savefig(F / "fig14_gc_survival.png", dpi=140); plt.close(fig)

# ---- fig15: forces on the strands ----
from ccsim.geometry import hopf_mirror_windings
from ccsim.forces import segment_forces
from ccsim.presets import PRESETS
pr = PRESETS["reactor"]()
w = hopf_mirror_windings(1.0, 0.005, 1.0, -1, eta=0.70, circuits=24, mode="woven", revolutions=2.0)
a = min(pr.wire_radius_m, 0.45 * w.clearance_m())
ws = w.with_current(pr.current_A).with_wire_radius(a)
fr = segment_forces(ws)
mid, fmag = fr["segments"]["mid"], np.linalg.norm(fr["segments"]["f"], axis=1)
fig = plt.figure(figsize=(11, 4.2))
ax0 = fig.add_subplot(1, 2, 1, projection="3d")
sc = ax0.scatter(mid[::3, 0], mid[::3, 1], mid[::3, 2], c=fmag[::3] / 1e6, s=1.5, cmap="magma")
ax0.set_title("|f| on the strands (MN/m) at 364 kA per strand", fontsize=9); ax0.set_box_aspect((1, 1, 0.6)); ax0.set_axis_off(); fig.colorbar(sc, ax=ax0, shrink=0.6, label="MN/m")
ax1 = fig.add_subplot(1, 2, 2)
ax1.hist(fmag / 1e6, bins=60, color="#d4884a")
ax1.set_xlabel("force per unit length (MN/m)"); ax1.set_ylabel("segments"); ax1.set_yscale("log"); ax1.grid(alpha=0.3)
ax1.set_title(f"peak {fr['f_max_N_m']/1e6:.1f} MN/m · B at conductor {fr['B_surface_max_T']:.1f} T · J {fr['J_A_m2']:.1e} A/m²", fontsize=9)
fig.tight_layout(); fig.savefig(F / "fig15_strand_forces.png", dpi=140); plt.close(fig)
json.dump({k: v for k, v in fr.items() if k != "segments"}, open(R / "exp12_forces.json", "w"), indent=1, default=float)
print("figures written")
