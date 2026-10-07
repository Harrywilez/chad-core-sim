"""Build agent/leaderboard.md from the ledger without re-evaluating anything: for every design (hash)
take its most recent evaluation with the default marker count.

    python3 agent/leaderboard_from_ledger.py [--since 2026-09-05]      # only rows evaluated after a date

Rows evaluated before the session-8 geometry correction (2026-09-05) describe a geometry that no
longer exists (touching strands, over-estimated clearance); use --since to exclude them.
"""

from __future__ import annotations

import argparse
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--since", default="2026-09-05")
    ap.add_argument("--particles", type=int, default=160)
    a = ap.parse_args()
    latest: dict = {}
    plasma: dict = {}
    for line in open(os.path.join(ROOT, "results", "design_ledger.jsonl")):
        r = json.loads(line)
        if r.get("time", "") < a.since or "particles" not in r:
            continue
        if r.get("n_particles", 160) != a.particles:
            continue
        # Fixed-current confinement studies have a different endpoint and no
        # legacy Score. Keep their observations out of this score leaderboard.
        if r.get("score") is None:
            continue
        h = r["hash"]
        if h not in latest or r["time"] > latest[h]["time"]:
            latest[h] = r
        if "score_plasma" in r and (h not in plasma or r["time"] > plasma[h]["time"]):
            plasma[h] = r
    rows = sorted(latest.values(), key=lambda r: -r["score"])
    lines = [f"# Leaderboard — score v2.2, corrected geometry (evaluations since {a.since}, {a.particles} markers, latest per design)", "",
             "score = τ_c × buildability, τ_c = √(τ_all·τ_pass) in wall transits at B_rms = 2 T; score_plasma applies the rung-3 polarization cap",
             "(v2.2: to every toroidal-field design without ≥ 20 % nested surfaces and |ι| above the electron-shorting threshold — closed-line",
             "designs and leaky toroidal solenoids alike; mirror designs, R ≥ 1.5, carry no factor). Every evaluation is in results/design_ledger.jsonl.", "",
             "| # | score | score_plasma | τ_all | τ_pass | build | clearance | S20 | S_end (T) | closed | R_mirror | design |", "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for i, r in enumerate(rows, 1):
        p, t, S = r["particles"], r["topology"], r["particles"]["survival"]
        s_end = S.get("S120", S.get("S80", S.get("S40", 0)))
        sp = plasma.get(r["hash"])
        sp_txt = f"{sp['score_plasma']:.1f}" + (" ✓" if sp["plasma"].get("surfaces_ok") else "") if sp else "–"
        lines.append(f"| {i} | {r['score']:.1f} | {sp_txt} | {p['tau_all']:.0f} | {p['tau_passing']:.0f} | {r['buildability']['factor']:.2f} | "
                     f"{r['validation']['clearance']:.4f} | {S.get('S20', 0):.2f} | {s_end:.2f} ({r.get('transits_run')}) | {t['frac_closed']:.2f} | "
                     f"{t['mirror_ratio_median']:.1f} | {r['name']}{' (invalid)' if r['invalid'] else ''} |")
    open(os.path.join(ROOT, "agent", "leaderboard.md"), "w").write("\n".join(lines) + "\n")
    print(f"{len(rows)} designs → agent/leaderboard.md")
    for l in lines[7:17]:
        print(l)


if __name__ == "__main__":
    main()
