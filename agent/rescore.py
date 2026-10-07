"""Re-score every design in designs/ with the current score and geometry, then write agent/leaderboard.md.

    python3 agent/rescore.py                 # all designs, 160 markers (≈ 30 s each)
    python3 agent/rescore.py --surfaces 3    # …and run --surfaces on the top 3 afterwards (≈ 2 min each)
    python3 agent/rescore.py --only mesh     # designs whose file name contains 'mesh'

Run it after any change to ccsim/ that moves numbers (geometry, score, seeds); the ledger keeps every
evaluation, so old and new rows can be compared by time.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
os.chdir(ROOT)
sys.path.insert(0, ROOT)

from ccsim.evaluate import evaluate  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--surfaces", type=int, default=0, help="run --surfaces (score_plasma) on the top N afterwards")
    ap.add_argument("--only", default="", help="substring filter on the design file name")
    ap.add_argument("--particles", type=int, default=160)
    a = ap.parse_args()
    rows = []
    files = sorted(f for f in glob.glob("designs/*.json") if a.only in os.path.basename(f))
    for f in files:
        d = json.load(open(f))
        try:
            r = evaluate(d, verbose=False, n_particles=a.particles)
        except Exception as e:  # noqa: BLE001
            print("ERR", f, e, flush=True)
            continue
        p = r.get("particles", {})
        S = p.get("survival", {})
        rows.append(dict(score=r["score"], tau_all=p.get("tau_all", 0), tau_pass=p.get("tau_passing", 0),
                         build=r.get("buildability", {}).get("factor", 0), S20=S.get("S20", 0),
                         S_end=S.get("S120", S.get("S80", S.get("S40", 0))), closed=r.get("topology", {}).get("frac_closed", 0),
                         T=r.get("transits_run"), file=os.path.basename(f), invalid=r["invalid"], clearance=r["validation"].get("clearance"),
                         design=d))
        print(f, round(r["score"], 1), "clearance", round(r["validation"].get("clearance", 0), 4), flush=True)
    rows.sort(key=lambda x: -x["score"])
    plasma = {}
    for row in rows[: a.surfaces]:
        try:
            r = evaluate(row["design"], verbose=False, surfaces=True, n_particles=a.particles)
            plasma[row["file"]] = (r.get("score_plasma"), r.get("surfaces", {}))
            print("surfaces", row["file"], r.get("score_plasma"), flush=True)
        except Exception as e:  # noqa: BLE001
            print("ERR surfaces", row["file"], e, flush=True)
    lines = ["# Leaderboard — score v2.1 (40/80/120-transit ladder, τ_c = √(τ_all·τ_pass), buildability on the conductor radius the model uses)", "",
             f"Re-scored with the current geometry (zone weave, exact clearance) and {a.particles} markers; every evaluation is also in results/design_ledger.jsonl.", ""]
    if plasma:
        lines += ["score_plasma (with --surfaces; electron-shorting criterion) for the finalists:", ""]
        for k, (sp, sf) in plasma.items():
            lines.append(f"* {k}: score_plasma {sp:.1f} — surfaces {100 * sf.get('frac_surface', 0):.0f} % of seeds, r_out {sf.get('outermost_surface_r', 0):.3f}, "
                         f"ι {sf.get('iota_axis', 0):.3f} → {sf.get('iota_edge', 0):.3f}, well {sf.get('well_depth', 0):+.2f}")
        lines.append("")
    lines += ["| # | score | τ_all | τ_pass | build | clearance | S20 | S_end | closed | run T | design |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for i, r in enumerate(rows, 1):
        lines.append(f"| {i} | {r['score']:.1f} | {r['tau_all']:.0f} | {r['tau_pass']:.0f} | {r['build']:.2f} | {r['clearance']:.4f} | {r['S20']:.2f} | {r['S_end']:.2f} | {r['closed']:.2f} | {r['T']} | {r['file']}{' (invalid)' if r['invalid'] else ''} |")
    open("agent/leaderboard.md", "w").write("\n".join(lines) + "\n")
    print("done:", len(rows), "designs")


if __name__ == "__main__":
    main()
