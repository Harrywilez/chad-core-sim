"""Figures for the playbook: fig09 (sweep ranking), fig10 (energy scan), fig11 (mirror-pair anatomy)."""

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
R, F = ROOT / "results", ROOT / "figures"
F.mkdir(exist_ok=True)

rows = json.load(open(R / "exp08_playbook.json"))
# ---- fig09: ranking ----
fig, ax = plt.subplots(figsize=(10, 12))
names = [r["name"] for r in rows][::-1]
y = np.arange(len(rows))
fam = lambda r: r["meta"].get("family", "")
colors = {"hopf-mirror": "#d4884a", "assembly": "#45b8a6", "hopf-torus": "#8a7bd6", "hopf-continued": "#7ad9ca", "hopf-drift": "#5b6873"}
for i, r in enumerate(rows[::-1]):
    c = colors.get(fam(r), "#e6b34f")
    ax.barh(i, r["S4"], color=c, alpha=0.25, height=0.8)
    ax.barh(i, r["S8"], color=c, alpha=0.5, height=0.8)
    ax.barh(i, r["S20"], color=c, alpha=1.0, height=0.8)
    if r["closed"] > 0.05:
        ax.text(r["S4"] + 0.01, i, f"closed lines {100*r['closed']:.0f}%", va="center", fontsize=7, color="#555")
ax.set_yticks(y); ax.set_yticklabels(names, fontsize=8)
ax.set_xlabel("reactor-scale D⁺ 15 keV survival: light = 4 transits, mid = 8, solid = 20")
ax.set_xlim(0, 1); ax.grid(axis="x", alpha=0.3)
ax.set_title("exp08 playbook — 42 configurations ranked by S20 (25³ grids, 160 particles, ±4 % statistical)")
from matplotlib.patches import Patch
ax.legend(handles=[Patch(color=v, label=k) for k, v in colors.items()] + [Patch(color="#e6b34f", label="new geometries")], loc="lower right", fontsize=8)
fig.tight_layout(); fig.savefig(F / "fig09_playbook_ranking.png", dpi=140); plt.close(fig)

# ---- fig10: energy scan ----
a = json.load(open(R / "exp09_a.json"))
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
series = {}
for key, rec in a.items():
    if rec["kind"] != "energy":
        continue
    series.setdefault(rec["name"], []).append((rec["run"]["energy_eV"] / 1e3, rec["run"]["rho_L_m"] * 1e3, rec["run"]["survival"]["S20"], rec["run"]["survival"]["S4"]))
for name, pts in series.items():
    pts.sort()
    E, rho, S20, S4 = zip(*pts)
    ax[0].plot(E, S20, "o-", label=name)
    ax[0].plot(E, S4, "o--", alpha=0.4, color=ax[0].lines[-1].get_color())
    ax[1].plot(rho, S20, "o-", label=name)
ax[0].set_xscale("log"); ax[0].set_xlabel("D⁺ energy (keV)"); ax[0].set_ylabel("survival (solid S20, dashed S4)"); ax[0].set_ylim(0, 0.8); ax[0].grid(alpha=0.3); ax[0].legend(fontsize=7)
ax[1].set_xlabel("Larmor radius at B_rms (mm, 1 m ball)"); ax[1].set_ylabel("S20"); ax[1].set_ylim(0, 0.8); ax[1].grid(alpha=0.3)
ax[1].set_title("drift-limited traps lose out as ρ_L grows; the meshed torus barely notices", fontsize=9)
fig.tight_layout(); fig.savefig(F / "fig10_energy_scan.png", dpi=140); plt.close(fig)

# ---- fig11: meshed torus anatomy (strands + |B| slices) ----
from ccsim.geometry import hopf_mirror_pair, hopf_mirror_windings
from ccsim.fields import biot_savart
A, B = hopf_mirror_pair(0.70, 24, 1.0, "woven", 0.025)
fig = plt.figure(figsize=(12, 4.2))
ax0 = fig.add_subplot(1, 3, 1, projection="3d")
ax0.plot(A[:, 0], A[:, 1], A[:, 2], lw=0.4, color="#d4884a"); ax0.plot(B[:, 0], B[:, 1], B[:, 2], lw=0.4, color="#4aa3d4")
ax0.set_title("torus strand (copper) + mirror strand (blue), woven", fontsize=9); ax0.set_box_aspect((1, 1, 0.6)); ax0.set_axis_off()
n = 121
xs = np.linspace(-0.87, 0.87, n)
for k, (sense, title) in enumerate(((-1, "opposed currents: toroidal field inside the tube"), (1, "same currents: poloidal (dipole-like) field"))):
    w = hopf_mirror_windings(1.0, 0.005, 1.0, sense, eta=0.70, circuits=24, mode="woven")
    X, Z = np.meshgrid(xs, xs, indexing="ij")
    P = np.column_stack((X.ravel(), np.zeros(X.size), Z.ravel()))
    Bf = biot_savart(P, w)
    Bm = np.linalg.norm(Bf, axis=1).reshape(n, n)
    axk = fig.add_subplot(1, 3, 2 + k)
    im = axk.imshow(np.log10(Bm.T + 1e-12), origin="lower", extent=(-0.87, 0.87, -0.87, 0.87), cmap="magma", vmin=-7.5, vmax=-4.5)
    axk.streamplot(xs, xs, Bf[:, 0].reshape(n, n).T, Bf[:, 2].reshape(n, n).T, color="w", linewidth=0.4, density=1.2, arrowsize=0.5)
    axk.add_patch(plt.Circle((0, 0), 0.82, fill=False, color="#e05d57", lw=0.8))
    axk.set_title(title, fontsize=9); axk.set_xlabel("x"); axk.set_ylabel("z")
    fig.colorbar(im, ax=axk, label="log10 |B| (T/A, unit ball)", shrink=0.8)
fig.tight_layout(); fig.savefig(F / "fig11_meshed_torus_anatomy.png", dpi=140); plt.close(fig)
print("figures written")
