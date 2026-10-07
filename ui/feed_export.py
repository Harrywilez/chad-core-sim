"""Export ledger rows as one JSON document each, for pushing into the published Lab's live
feed (the artifact's shared database, collection `designs`, doc id = row id).

    python3 ui/feed_export.py                # every row → cache/feed/<id>.json (skips ones already exported)
    python3 ui/feed_export.py --after 240    # only ledger lines after 240
    python3 ui/feed_export.py --status "run 3: tornado plugs" --running 1 --note "…"   # status/run document

The artifact database is written by Claude (Artifact tool, write_db batch with file_path per
document), so this is the hand-off format; the local live view (ui/serve.py) needs none of it.
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "ui"))
from ledger_rows import read_rows  # noqa: E402

OUT = ROOT / "cache" / "feed"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--after", type=int, default=0)
    ap.add_argument("--status", default=None, help="label for the status/run document")
    ap.add_argument("--running", type=int, default=None)
    ap.add_argument("--note", default="")
    ap.add_argument("--best", default="")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    rows, n = read_rows(a.after)
    written = []
    for r in rows:
        f = OUT / f"{r['id']}.json"
        if f.exists() and a.after == 0:
            continue
        json.dump(r, open(f, "w"), default=float)
        written.append(f)
    print(f"ledger lines: {n}; exported {len(written)} new rows to {OUT}")
    for f in written:
        print("  ", f.name)
    if a.status is not None:
        st = {"label": a.status, "running": bool(a.running), "note": a.note, "best": a.best, "evaluations": n,
              "updated": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")}
        json.dump(st, open(OUT / "status_run.json", "w"))
        print("status:", st)


if __name__ == "__main__":
    main()
