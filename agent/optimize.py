"""Component composition search with parameter tuning on an existing design JSON.

    python3 agent/optimize.py designs/seed.json --iters 30 [--sigma 0.15] [--out designs/seed_opt.json]

Mutations (one or two per step): scale a numeric parameter by (1 ± σ·N(0,1)), nudge a
translation by σ·0.3, rotate by σ·40°, flip a component current, change the wire radius by
±20 %.  A candidate is accepted when its score improves; every evaluation is appended to
results/design_ledger.jsonl by ccsim.evaluate.  Structural add/remove/replace/duplicate moves are enabled by default; --parameter-only
retains the original hill climb. --library supplies donor components from design JSONs.
Valid lower-scoring candidates enter a bounded exploration population.
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import uuid
import random
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from agent.composition import propose, component_library, choose_parent  # noqa: E402

INT_PARAMS = {"circuits", "turns", "periods", "rings"}
FROZEN = {"variant", "mode", "helical_mode", "axis", "alternate", "alternating", "sense", "levels"}


def mutate(design: dict, sigma: float, rng: random.Random) -> dict:
    d = copy.deepcopy(design)
    comps = d["components"]
    for _ in range(rng.choice((1, 1, 2))):
        c = rng.choice(comps)
        move = rng.random()
        if move < 0.45 and c.get("params"):
            keys = [k for k, v in c["params"].items() if isinstance(v, (int, float)) and not isinstance(v, bool) and k not in FROZEN]
            if keys:
                k = rng.choice(keys)
                v = c["params"][k]
                nv = v * (1 + sigma * rng.gauss(0, 1))
                c["params"][k] = max(1, int(round(nv))) if k in INT_PARAMS else round(nv, 4)
        elif move < 0.65:
            c["scale"] = round(max(0.1, c.get("scale", 1.0) * (1 + sigma * rng.gauss(0, 1))), 4)
        elif move < 0.8:
            t = list(c.get("translate", [0.0, 0.0, 0.0]))
            j = rng.randrange(3)
            t[j] = round(t[j] + 0.3 * sigma * rng.gauss(0, 1), 4)
            c["translate"] = t
        elif move < 0.9:
            r = c.get("rotate", {"axis": [0, 0, 1], "deg": 0.0})
            c["rotate"] = {"axis": r.get("axis", [0, 0, 1]), "deg": round(r.get("deg", 0.0) + 40 * sigma * rng.gauss(0, 1), 2)}
        elif move < 0.95:
            c["current"] = -c.get("current", 1.0)
        else:
            d["wire_radius"] = round(max(0.0003, d.get("wire_radius", 0.005) * (1.2 if rng.random() < 0.5 else 0.8)), 5)
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("seed_design")
    ap.add_argument("--iters", type=int, default=30)
    ap.add_argument("--sigma", type=float, default=0.15)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out")
    ap.add_argument("--parameter-only", action="store_true", help="legacy fixed-structure hill climb")
    ap.add_argument("--library", nargs="*", default=[], help="donor design JSON files")
    ap.add_argument("--max-components", type=int, default=6)
    ap.add_argument("--structural-rate", type=float, default=0.5)
    ap.add_argument("--exploration", type=float, default=0.25)
    ap.add_argument("--population", type=int, default=24)
    ap.add_argument("--position-seed", type=int, default=271828)
    ap.add_argument("--topology-seed", type=int, default=31415)
    ap.add_argument("--history", default="results/composition_search.jsonl")
    ap.add_argument("--particles", type=int, default=160)
    ap.add_argument("--grid", type=int, default=25)
    ap.add_argument("--transits", type=int, default=40)
    ap.add_argument("--max-transits", type=int, default=120)
    ap.add_argument("--eval-seed", type=int, default=12345)
    ap.add_argument("--surfaces", action="store_true", help="optimize score_plasma; slower")
    a = ap.parse_args()
    from ccsim.evaluate import evaluate
    from ccsim.design import FAMILIES
    if a.iters < 0 or a.max_components < 1 or not 0 <= a.structural_rate <= 1 or not 0 <= a.exploration <= 1 or a.population < 1:
        ap.error("invalid iteration count, component cap, or structural rate")
    if a.particles < 1 or a.grid < 3 or a.transits < 1 or a.max_transits < max(40, a.transits) or a.sigma < 0:
        ap.error("particles/grid/sigma invalid or max-transits below the evaluator minimum of 40")
    rng = random.Random(a.seed)
    best = json.loads(Path(a.seed_design).read_text())
    base_name = best.get("name", "design")
    if not best.get("components") or (not a.parameter_only and len(best["components"]) > a.max_components):
        ap.error("seed component count must be in [1,max-components]")
    donors = component_library([json.loads(Path(p).read_text()) for p in a.library])
    if not donors:
        donors = [{"family": name, "params": copy.deepcopy(info["params"]),
                   "scale": 0.35, "translate": [0., 0., 0.45], "current": 0.3}
                  for name, info in FAMILIES.items()]
    options = dict(n_particles=a.particles, grid_n=a.grid, transits=a.transits,
                   max_transits=a.max_transits, seed=a.eval_seed, surfaces=a.surfaces,
                   position_seed=a.position_seed, topology_seed=a.topology_seed)
    metric = "score_plasma" if a.surfaces else "score"
    run_id = uuid.uuid4().hex
    if a.out:
        Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    history = Path(a.history)
    history.parent.mkdir(parents=True, exist_ok=True)
    def assess(candidate, iteration, move):
        try:
            result = evaluate(candidate, verbose=True, **options)
        except (ValueError, KeyError, OverflowError, ZeroDivisionError) as exc:
            result = {"invalid": True, "error": repr(exc), "score": 0.0}
        value = result.get(metric)
        # A failed surface analysis cannot silently compete on raw score.
        score = float(value) if value is not None and not result.get("invalid") else float("-inf")
        if not math.isfinite(score):
            score = float("-inf")
        record = dict(run_id=run_id, iteration=iteration, move=move, search_seed=a.seed,
                      search_options=vars(a), donors=donors if iteration == -1 else None,
                      evaluation_options=options, metric=metric, design=candidate,
                      result=result, objective=score if score != float("-inf") else None)
        with history.open("a") as f:
            f.write(json.dumps(record) + "\n")
        return score
    best_score = assess(best, -1, "baseline")
    print(f"start: score {best_score:.2f}")
    population = [copy.deepcopy(best)] if best_score != float("-inf") else []
    for it in range(a.iters):
        if a.parameter_only:
            cand, move = mutate(best, a.sigma, rng), "parameter"
        else:
            parent = choose_parent(best, population, rng, a.exploration)
            cand, move = propose(parent, rng, donors, mutate, a.sigma,
                                 a.max_components, a.structural_rate)
        cand["name"] = f"{base_name} · opt{it}"
        score = assess(cand, it, move)
        if score != float("-inf"):
            population.append(copy.deepcopy(cand))
            population = population[-a.population:]
        if score > best_score:
            best, best_score = cand, score
            print(f"  ↑ accepted at iteration {it}: score {best_score:.2f}")
            if a.out:
                Path(a.out).write_text(json.dumps(best, indent=1) + "\n")
    if not math.isfinite(best_score):
        raise SystemExit("No valid candidate with a finite requested objective; inspect history")
    print(f"best {metric} {best_score:.2f}")
    if a.out:
        Path(a.out).write_text(json.dumps(best, indent=1) + "\n")
    print(json.dumps(best, indent=1))


if __name__ == "__main__":
    main()
