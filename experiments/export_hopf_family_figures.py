"""Figures for exports/hopf_family/README.md (run after export_hopf_family.py).

    python3 experiments/export_hopf_family_figures.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from ccsim import geometry as G  # noqa: E402

OUT = ROOT / "exports" / "hopf_family"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)
M = json.load(open(OUT / "MANIFEST.json"))
BLUE, RED, GREEN = "#2a73d9", "#d94d33", "#33905a"
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "figure.dpi": 130})


def pts(path):
    """x,y,z CSV with '#' comment lines and an x,y,z header."""
    rows = [l for l in open(path, encoding="utf-8") if l.strip() and not l.startswith("#")]
    assert rows[0].strip() == "x,y,z", rows[0]
    return np.array([[float(v) for v in l.split(",")] for l in rows[1:]])


def torus_surface(R0, r, n=60, m=30):
    th, ph = np.meshgrid(np.linspace(0, 2 * np.pi, m), np.linspace(0, 2 * np.pi, n), indexing="ij")
    return (R0 + r * np.cos(th)) * np.cos(ph), (R0 + r * np.cos(th)) * np.sin(ph), r * np.sin(th)


def set_equal(ax, lim=1.0):
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_zlim(-lim, lim)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xticks([-1, 0, 1]); ax.set_yticks([-1, 0, 1]); ax.set_zticks([-1, 0, 1])


# ---------------------------------------------------------------------------
# Fig 1: Hopf coordinates — base sphere paths and Villarceau circles
# ---------------------------------------------------------------------------


def fig1():
    fig = plt.figure(figsize=(11, 5))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    # base sphere (radius 1/2): point (η, φ) ↦ polar angle 2η, azimuth φ
    th, ph = np.meshgrid(np.linspace(0, np.pi, 40), np.linspace(0, 2 * np.pi, 80), indexing="ij")
    ax.plot_surface(0.5 * np.sin(th) * np.cos(ph), 0.5 * np.sin(th) * np.sin(ph), 0.5 * np.cos(th), color="#dddddd", alpha=0.25, linewidth=0)
    for eta, col, lab in ((0.82, "#888888", "η = 0.82 (start torus)"), (0.70, "#555555", "η = 0.70 (torus of v1)"), (0.60, "#222222", "η = 0.60 (torus of II)")):
        p = np.linspace(0, 2 * np.pi, 200)
        ax.plot(0.5 * np.sin(2 * eta) * np.cos(p), 0.5 * np.sin(2 * eta) * np.sin(p), np.full_like(p, 0.5 * np.cos(2 * eta)), color=col, lw=1.2, label=lab)
    cols = plt.cm.viridis(np.linspace(0.1, 0.9, 6))
    for (key, spec), col in zip(G.HOPF_DRIFT_VARIANTS.items(), cols):
        v = np.linspace(0, 1, 100)
        eta = 0.82 + (spec["eta1"] - 0.82) * v
        phi = math.radians(spec["sweep"]) * v
        ax.plot(0.5 * np.sin(2 * eta) * np.cos(phi), 0.5 * np.sin(2 * eta) * np.sin(phi), 0.5 * np.cos(2 * eta), color=col, lw=2, label=f"drift {key}")
    ax.scatter([0.5 * np.sin(1.64)], [0], [0.5 * np.cos(1.64)], color="k", s=20)
    ax.set_title("Hopf base sphere: a torus is a latitude circle,\nthe drift winding's fibre walks across it")
    ax.set_xlim(-0.55, 0.55); ax.set_ylim(-0.55, 0.55); ax.set_zlim(-0.55, 0.55); ax.set_box_aspect((1, 1, 1))
    ax.set_xticks([]); ax.set_yticks([]); ax.set_zticks([])
    ax.legend(loc="upper left", fontsize=7, bbox_to_anchor=(-0.15, 1.0))
    ax.view_init(elev=22, azim=-50)

    ax = fig.add_subplot(1, 2, 2, projection="3d")
    eta = 0.70
    fr = M["facts"]["torus_frames"]["0.70"]
    s = 1.0 / (fr["R0_raw"] + fr["r_raw"])
    X, Y, Z = torus_surface(fr["R0"], fr["r"])
    ax.plot_surface(X, Y, Z, color="#cccccc", alpha=0.18, linewidth=0)
    u = np.linspace(0, 1, 400)
    for k in range(4):
        phi = 2 * np.pi * k / 4
        f = G._hopf_points(np.full_like(u, eta), np.full_like(u, phi), 2 * np.pi * u) * s
        ax.plot(f[:, 0], f[:, 1], f[:, 2], color=BLUE, lw=1.4, label="Hopf fibres (family A)" if k == 0 else None)
        g = f * np.array([1, 1, -1.0])
        ax.plot(g[:, 0], g[:, 1], g[:, 2], color=RED, lw=1.4, label="their z-mirror images (family B)" if k == 0 else None)
    ax.set_title("Torus η = 0.70: fibres are Villarceau circles;\nthe mirror image of a fibre is an anti-fibre — each pair meets twice")
    set_equal(ax, 1.05)
    ax.legend(loc="upper left", fontsize=7)
    ax.view_init(elev=32, azim=-55)
    fig.tight_layout()
    fig.savefig(FIG / "fig1_hopf_coordinates.png")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Fig 2: the six drift variants + the continuation
# ---------------------------------------------------------------------------


def fig2():
    fig = plt.figure(figsize=(12, 8))
    keys = list(G.HOPF_DRIFT_VARIANTS)
    for i, key in enumerate(keys):
        ax = fig.add_subplot(2, 3, i + 1, projection="3d")
        P = pts(OUT / "1_hopf_drift" / f"hopf_drift_{key}_10circuits.drift.csv")
        num = M["geometries"][f"hopf_drift_{key}"]
        eta, a1, a2 = G.hopf_coordinates(P, num["scale_s0"])
        seg = np.stack((P[:-1], P[1:]), axis=1)
        from mpl_toolkits.mplot3d.art3d import Line3DCollection
        lc = Line3DCollection(seg, cmap="viridis", norm=plt.Normalize(0.48, 0.82), linewidths=1.4)
        lc.set_array(0.5 * (eta[:-1] + eta[1:]))
        ax.add_collection(lc)
        ax.scatter(*P[0], color="k", s=14); ax.scatter(*P[-1], color="k", marker="x", s=20)
        set_equal(ax, 1.0)
        ax.set_title(f"{key}: η 0.82 → {num['eta1']}, sweep {num['sweep_deg']:.0f}°\nL = {num['length']:.1f}, clearance {num['self_clearance']:.3f}")
        ax.view_init(elev=30, azim=-60)
    fig.suptitle("Hopf-drift windings (10 fibre circuits, colour = η: yellow fat torus → purple thin torus; ● start, × end)")
    fig.tight_layout()
    fig.savefig(FIG / "fig2_hopf_drift_variants.png")
    plt.close(fig)

    fig = plt.figure(figsize=(12, 4))
    P = pts(OUT / "1_hopf_drift" / "hopf_drift_s64_180_10circuits.drift.csv")
    for j, (el, az, ttl) in enumerate(((90, -90, "s64_180 — top view"), (0, -90, "s64_180 — side view"))):
        ax = fig.add_subplot(1, 3, j + 1, projection="3d")
        ax.plot(P[:, 0], P[:, 1], P[:, 2], color=GREEN, lw=1.0)
        ax.scatter(*P[0], color="k", s=14); ax.scatter(*P[-1], color="k", marker="x", s=20)
        set_equal(ax, 1.0); ax.view_init(elev=el, azim=az); ax.set_title(ttl)
    ax = fig.add_subplot(1, 3, 3, projection="3d")
    Pc = pts(OUT / "1_hopf_drift" / "hopf_continued_720deg_48circuits.drift.csv")
    ax.plot(Pc[:, 0], Pc[:, 1], Pc[:, 2], color=GREEN, lw=0.5)
    set_equal(ax, 1.0); ax.view_init(elev=30, azim=-60)
    ax.set_title("continuation: φ sweeps 720°, 48 circuits —\nthe drift keeps going and tiles the torus")
    fig.tight_layout()
    fig.savefig(FIG / "fig2b_hopf_drift_views.png")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Helpers for the pair figures
# ---------------------------------------------------------------------------


def poloidal_section(P, phi0, R0, r, width=0.02):
    """Points of a polyline within a thin wedge |φ − φ0| < width, mapped to (ρ − R0, z) of the poloidal plane."""
    phi = np.arctan2(P[:, 1], P[:, 0])
    d = np.angle(np.exp(1j * (phi - phi0)))
    sel = np.abs(d) < width
    rho = np.hypot(P[sel, 0], P[sel, 1])
    return rho - R0, P[sel, 2]


def unrolled(P, scale, a1_win, a2_win):
    """(a1, a2) of the points inside a window, plus recovered η, with line breaks at the wrap."""
    eta, a1, a2 = G.hopf_coordinates(P, scale)
    a1 = np.mod(a1, 2 * np.pi); a2 = np.mod(a2, 2 * np.pi)
    return eta, a1, a2


def draw_unrolled(ax, A, B, scale, eta0, de, win=(0.0, 1.2, 0.0, 1.2), title=""):
    eA, a1A, a2A = unrolled(A, scale, None, None)
    eB, a1B, a2B = unrolled(B, scale, None, None)
    for (a1, a2, e, col, lab) in ((a1A, a2A, eA, BLUE, "A"), (a1B, a2B, eB, RED, "B")):
        # break the polyline at wraps
        d = np.hypot(np.diff(a1), np.diff(a2))
        brk = np.where(d > 1.0)[0]
        segs = np.split(np.arange(len(a1)), brk + 1)
        for s in segs:
            if len(s) < 2:
                continue
            x, y = a1[s], a2[s]
            inwin = (x >= win[0] - 0.05) & (x <= win[1] + 0.05) & (y >= win[2] - 0.05) & (y <= win[3] + 0.05)
            if not inwin.any():
                continue
            seg = np.stack((np.column_stack((x[:-1], y[:-1])), np.column_stack((x[1:], y[1:]))), axis=1)
            lw = 0.9 + 1.6 * (0.5 * (e[s][:-1] + e[s][1:]) - eta0) / de      # thicker = displaced outward (over)
            lc = LineCollection(seg, colors=col, linewidths=np.clip(lw, 0.3, 3.0), alpha=0.95)
            ax.add_collection(lc)
    ax.set_xlim(win[0], win[1]); ax.set_ylim(win[2], win[3])
    ax.set_xlabel("a₁ (toroidal angle, rad)"); ax.set_ylabel("a₂ (poloidal, rad)")
    ax.set_title(title)
    ax.set_aspect("equal")


def draw_weave_patch(ax, A, B, scale, eta0, de, centre_a1, centre_a2, half=0.35, title="", lw_base=5.0):
    """A patch of the torus seen from outside along the surface normal, drawn in (a₁, a₂) with proper
    occlusion: strands are sorted by their recovered η (outer drawn last) so over/under is shown as it is."""
    eA, a1A, a2A = G.hopf_coordinates(A, scale)
    eB, a1B, a2B = G.hopf_coordinates(B, scale)
    pieces = []
    for (a1, a2, e, col) in ((a1A, a2A, eA, BLUE), (a1B, a2B, eB, RED)):
        x = np.angle(np.exp(1j * (a1 - centre_a1)))
        y = np.angle(np.exp(1j * (a2 - centre_a2)))
        inside = (np.abs(x) < half * 1.15) & (np.abs(y) < half * 1.15)
        idx = np.where(inside)[0]
        runs = np.split(idx, np.where(np.diff(idx) > 1)[0] + 1)         # contiguous stretches inside the patch
        for rr in runs:
            if len(rr) < 2:
                continue
            side = np.sign(e[rr] - eta0)
            cut = np.where(np.diff(side) != 0)[0] + 1                    # split where the strand goes from over to under
            for pc in np.split(rr, cut):
                if len(pc) >= 1:
                    pc2 = np.r_[pc, min(pc[-1] + 1, len(e) - 1)]           # overlap one sample so pieces join
                    pieces.append((float(e[pc].mean()), x[pc2], y[pc2], col))
    pieces.sort(key=lambda t: t[0])                                       # draw the innermost first, outermost last
    for e, x, y, col in pieces:
        ax.plot(x, y, color="white", lw=lw_base + 2.4, solid_capstyle="butt", solid_joinstyle="round")
        ax.plot(x, y, color=col, lw=lw_base, solid_capstyle="butt", solid_joinstyle="round")
    ax.set_xlim(-half, half); ax.set_ylim(-half, half); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title)


def fig_pair(key_csv_dir, key, scale, eta0, de, title, fname, closed, extra_panel=None):
    A = pts(key_csv_dir / f"{key}.A.csv")
    B = pts(key_csv_dir / f"{key}.B.csv")
    fr = M["facts"]["torus_frames"][f"{eta0:.2f}"]
    R0, r = fr["R0"], fr["r"]
    fig = plt.figure(figsize=(13, 8.5))
    ax = fig.add_subplot(2, 3, 1, projection="3d")
    ax.plot(A[:, 0], A[:, 1], A[:, 2], color=BLUE, lw=0.5, label="strand A")
    ax.plot(B[:, 0], B[:, 1], B[:, 2], color=RED, lw=0.5, label="strand B")
    if not closed:
        ax.scatter(*A[0], color=BLUE, s=25, edgecolor="k"); ax.scatter(*A[-1], color=BLUE, marker="x", s=30)
        ax.scatter(*B[0], color=RED, s=25, edgecolor="k"); ax.scatter(*B[-1], color=RED, marker="x", s=30)
    set_equal(ax, 1.0); ax.view_init(elev=32, azim=-55); ax.legend(loc="upper left", fontsize=7)
    ax.set_title("perspective" + ("" if closed else " (● start, × end)"))
    ax = fig.add_subplot(2, 3, 2, projection="3d")
    ax.plot(A[:, 0], A[:, 1], A[:, 2], color=BLUE, lw=0.4)
    ax.plot(B[:, 0], B[:, 1], B[:, 2], color=RED, lw=0.4)
    set_equal(ax, 1.0); ax.view_init(elev=90, azim=-90); ax.set_title("top view")
    # poloidal section
    ax = fig.add_subplot(2, 3, 3)
    for phi0 in (0.3,):
        xa, za = poloidal_section(A, phi0, R0, r)
        xb, zb = poloidal_section(B, phi0, R0, r)
        ax.scatter(xa, za, s=9, color=BLUE, label="A")
        ax.scatter(xb, zb, s=9, color=RED, label="B")
    th = np.linspace(0, 2 * np.pi, 300)
    ax.plot(r * np.cos(th), r * np.sin(th), color="#999999", lw=0.8, label=f"torus η = {eta0}")
    ax.set_aspect("equal"); ax.set_xlabel("ρ − R₀"); ax.set_ylabel("z"); ax.legend(fontsize=7, loc="upper right")
    ax.set_title("poloidal cross-section (points where the strands\npierce the plane φ = 0.3 ± 0.02)")
    # unrolled window
    ax = fig.add_subplot(2, 3, 4)
    draw_unrolled(ax, A, B, scale, eta0, de, win=(0.0, 1.0, 0.0, 1.0), title="unrolled (a₁, a₂): A slope +q/p, B slope −q/p;\nline width = radial displacement (thick = out)")
    # recovered η along A
    ax = fig.add_subplot(2, 3, 5)
    eA, a1A, a2A = G.hopf_coordinates(A, scale)
    eB, a1B, a2B = G.hopf_coordinates(B, scale)
    s = np.linspace(0, 1, len(A))
    n = int(0.025 * len(A))
    ax.plot(s[:n], eA[:n], color=BLUE, lw=1, label="A")
    ax.plot(s[:n], eB[:n], color=RED, lw=1, label="B")
    ax.axhline(eta0, color="#999999", lw=0.6); ax.axhline(eta0 + de, color="#cccccc", lw=0.6, ls="--"); ax.axhline(eta0 - de, color="#cccccc", lw=0.6, ls="--")
    ax.set_xlabel("fraction of strand length"); ax.set_ylabel("recovered η"); ax.legend(fontsize=7)
    ax.set_title("radial weave profile η(s) along the first 2.5 % of each strand\n(η ± δη = over/under; ramps between crossings)")
    ax = fig.add_subplot(2, 3, 6)
    if extra_panel is not None:
        extra_panel(ax, A, B)
    else:
        draw_weave_patch(ax, A, B, scale, eta0, de, centre_a1=0.9, centre_a2=1.5708, half=0.32, title="weave close-up on the outboard side (seen from outside;\nblue A, red B, drawn in over/under order)")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(FIG / fname)
    plt.close(fig)


def fig3():
    m = M["geometries"]["mirror_pair_v1_woven"]
    fig_pair(OUT / "2_mirror_pair_v1", "mirror_pair_v1_woven", m["strand_scale"], 0.70, 0.025,
             f"Mirror pair v1 (session-5 scheme): η 0.70, N 24, 1 revolution, δη 0.025 — open (24½, 23½) strands, B = z-mirror, {m['crossings']} crossings, min distance {m['min_distance_A_B']:.4f}",
             "fig3_mirror_pair_v1.png", closed=False)


def fig4():
    m = M["geometries"]["meshed_torus_II_achiral_de015"]
    fig_pair(OUT / "3_meshed_torus_II", "meshed_torus_II_achiral_de015", m["strand_scale"], 0.60, 0.015,
             f"Meshed torus II, achiral: η 0.60, N 24, 2 revolutions, δη 0.015 — closed (25, 23) knots, B = −A, {m['crossings']} crossings, zone weave ({m['weave_regime']}, {100*m['zone_fraction']:.0f} % frozen zones), min distance {m['min_distance_A_B']:.4f}",
             "fig4_meshed_torus_II_achiral.png", closed=True)


def fig5():
    m = M["geometries"]["meshed_torus_II_chiral_e070_n4_de015"]
    key = "meshed_torus_II_chiral_e070_n4_de015"
    A = pts(OUT / "3_meshed_torus_II" / f"{key}.A.csv")
    B = pts(OUT / "3_meshed_torus_II" / f"{key}.B.csv")
    scale, eta0, de = m["strand_scale"], 0.60, 0.015
    fr = M["facts"]["torus_frames"]["0.60"]
    R0, r = fr["R0"], fr["r"]
    fig = plt.figure(figsize=(13, 8.5))
    ax = fig.add_subplot(2, 3, 1, projection="3d")
    ax.plot(A[:, 0], A[:, 1], A[:, 2], color=BLUE, lw=0.5, label="A")
    ax.plot(B[:, 0], B[:, 1], B[:, 2], color=RED, lw=0.5, label="B")
    set_equal(ax, 1.0); ax.view_init(elev=32, azim=-55); ax.legend(loc="upper left", fontsize=7); ax.set_title("perspective")
    ax = fig.add_subplot(2, 3, 2, projection="3d")
    ax.plot(A[:, 0], A[:, 1], A[:, 2], color=BLUE, lw=0.4)
    ax.plot(B[:, 0], B[:, 1], B[:, 2], color=RED, lw=0.4)
    set_equal(ax, 1.0); ax.view_init(elev=90, azim=-90); ax.set_title("top view: the n = 4 modulation shows as four-fold bunching")
    # unrolled: turn density
    ax = fig.add_subplot(2, 3, 3)
    draw_unrolled(ax, A, B, scale, eta0, de, win=(0.0, 3.2, 0.0, 6.3), title="unrolled: turns bunch where 1 + ε cos(2a₂ − 4a₁) is large\n(both families the same-handed → chiral)")
    # theoretical density overlay
    a1g, a2g = np.meshgrid(np.linspace(0, 3.2, 200), np.linspace(0, 6.3, 300))
    ax.contour(a1g, a2g, 1 + 0.7 * np.cos(2 * a2g - 4 * a1g), levels=[1.0], colors="#444444", linewidths=0.8, linestyles="--")
    # two cross-sections half a period apart
    for j, phi0 in enumerate((0.0, np.pi / 4)):
        ax = fig.add_subplot(2, 3, 4 + j)
        xa, za = poloidal_section(A, phi0, R0, r)
        xb, zb = poloidal_section(B, phi0, R0, r)
        ax.scatter(xa, za, s=9, color=BLUE); ax.scatter(xb, zb, s=9, color=RED)
        th = np.linspace(0, 2 * np.pi, 300)
        ax.plot(r * np.cos(th), r * np.sin(th), color="#999999", lw=0.8)
        # density arrows: where is the conductor dense? mark the poloidal angle of max density 2a2 - 4 phi = 0
        ax.set_aspect("equal"); ax.set_xlabel("ρ − R₀"); ax.set_ylabel("z")
        ax.set_title(f"cross-section at φ = {phi0:.2f} rad ({'0' if j == 0 else 'π/4 = half a period'}):\nthe dense side rotates with φ (l = 2 pattern)")
    ax = fig.add_subplot(2, 3, 6)
    draw_weave_patch(ax, A, B, scale, eta0, de, centre_a1=0.9, centre_a2=1.5708, half=0.32, title="weave close-up (outboard side)")
    fig.suptitle(f"Meshed torus II, chiral: η 0.60, N 24, ε 0.70, n 4, δη 0.015 — {m['crossings']} crossings, {m['weave_regime']} weave ({100*m['zone_fraction']:.0f} % zones), min distance {m['min_distance_A_B']:.4f}; the design with nested flux surfaces")
    fig.tight_layout()
    fig.savefig(FIG / "fig5_meshed_torus_II_chiral.png")
    plt.close(fig)


def fig6():
    fig, axs = plt.subplots(1, 4, figsize=(14, 3.9))
    items = [("2_mirror_pair_v1", "mirror_pair_v1_woven", 0.70, 0.025, "v1: δη 0.025, lattice weave (−1)^k"),
             ("3_meshed_torus_II", "meshed_torus_II_achiral_de015", 0.60, 0.015, "II achiral: δη 0.015, zone weave (mixed)"),
             ("3_meshed_torus_II", "meshed_torus_II_achiral_de010", 0.60, 0.010, "II achiral: δη 0.010, zone weave (woven)"),
             ("3_meshed_torus_II", "meshed_torus_II_chiral_e070_n4_de015", 0.60, 0.015, "II chiral ε 0.7 n 4: δη 0.015")]
    for ax, (d, key, eta0, de, ttl) in zip(axs, items):
        A = pts(OUT / d / f"{key}.A.csv"); B = pts(OUT / d / f"{key}.B.csv")
        draw_weave_patch(ax, A, B, M["geometries"][key]["strand_scale"], eta0, de, centre_a1=0.9, centre_a2=1.5708, half=0.30, title=ttl, lw_base=4.5)
    fig.suptitle("Weave close-ups on the outboard side, same patch of (a₁, a₂) for each; strands drawn in over/under order")
    fig.tight_layout()
    fig.savefig(FIG / "fig6_weave_closeups.png")
    plt.close(fig)


if __name__ == "__main__":
    for f in (fig1, fig2, fig3, fig4, fig5, fig6):
        f()
        print("done", f.__name__)
