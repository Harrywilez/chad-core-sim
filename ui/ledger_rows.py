"""Compact rows for the Lab's Designs panel, from `results/design_ledger.jsonl`.

One evaluation → one small JSON object (the design spec plus the numbers the panel shows).
Used by build.py (embedded snapshot), serve.py (local live feed) and feed_export.py (artifact feed).
"""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LEDGER = ROOT / "results" / "design_ledger.jsonl"


def family_summary(design: dict) -> str:
    fams: dict = {}
    for c in design.get("components", []):
        fams[c.get("family", "?")] = fams.get(c.get("family", "?"), 0) + 1
    return " + ".join(f"{f}×{n}" if n > 1 else f for f, n in fams.items())


def compact(r: dict, line_no: int | None = None) -> dict:
    p = r.get("particles") or {}
    S = p.get("survival") or {}
    S_end = None
    for k in ("S120", "S80", "S40", "S20"):
        if k in S:
            S_end = S[k]
            break
    f = r.get("field") or {}
    b = r.get("buildability") or {}
    t = r.get("topology") or {}
    v = r.get("validation") or {}
    row = {
        "id": f"{r.get('hash', 'x')}_{(r.get('time') or '').replace(':', '').replace('-', '')}",
        "hash": r.get("hash"),
        "name": r.get("name"),
        "t": r.get("time"),
        "fam": family_summary(r.get("design") or {}),
        "score": r.get("score"),
        "tau_c": r.get("tau_c"),
        "tau_all": p.get("tau_all"),
        "tau_pass": p.get("tau_passing"),
        "build": b.get("factor"),
        "closed": t.get("frac_closed"),
        "pred": t.get("predicted_adiabatic_retention"),
        "S4": S.get("S4"),
        "S20": S.get("S20"),
        "S_end": S_end,
        "T": r.get("transits_run") or r.get("transits"),
        "n": r.get("n_particles"),
        "I_kA": (f.get("current_A") or 0) / 1e3 if f else None,
        "brms364": f.get("B_rms_per_364kA_T"),
        "clr": v.get("clearance"),
        "wire": b.get("wire_radius_used"),
        "valid": not r.get("invalid", False),
        "problems": v.get("problems") or [],
        "design": r.get("design"),
        "secs": r.get("seconds"),
    }
    if line_no is not None:
        row["line"] = line_no
    if r.get("surfaces"):
        s = r["surfaces"]
        row["surf"] = {"surf": s.get("frac_surface"), "r_out": s.get("outermost_surface_r"), "iota_axis": s.get("iota_axis"),
                       "iota_edge": s.get("iota_edge"), "well": s.get("well_depth")}
        row["score_plasma"] = r.get("score_plasma")
        pl = (r.get("plasma") or {}).get("polarization") or {}
        row["plasma"] = {"tau_pol": (r.get("plasma") or {}).get("tau_pol_wall_transits"), "ok": (r.get("plasma") or {}).get("surfaces_ok"),
                         "factor": (r.get("plasma") or {}).get("factor"), "iota_needed": pl.get("iota_needed_to_short")}
    return _clean(row)


def _clean(x):
    """NaN/inf → null (strict JSON for the feed and the page)."""
    if isinstance(x, float):
        return x if x == x and abs(x) != float("inf") else None
    if isinstance(x, dict):
        return {k: _clean(v) for k, v in x.items()}
    if isinstance(x, list):
        return [_clean(v) for v in x]
    return x


def read_rows(after: int = 0, path: Path = LEDGER) -> tuple[list, int]:
    """Compact rows for ledger lines after `after` (0-based count), and the new line count."""
    if not path.exists():
        return [], 0
    out = []
    n = 0
    with open(path, "r", encoding="utf-8") as fh:
        for n, line in enumerate(fh, start=1):
            if n <= after:
                continue
            line = line.strip()
            if not line:
                continue
            try:
                out.append(compact(json.loads(line), n))
            except Exception as e:  # a half-written last line while an evaluation is appending
                if n == after + len(out) + 1:
                    return out, n - 1
                print("ledger line", n, "skipped:", e)
    return out, n


if __name__ == "__main__":
    rows, n = read_rows()
    print(n, "lines;", len(rows), "rows; last:", json.dumps(rows[-1], default=float)[:300] if rows else None)
