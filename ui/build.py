"""Assemble the Chad Core Lab page from its parts.

  chad_core_lab.html            artifact body (no doctype/html/head/body — the host adds them)
  chad_core_lab_standalone.html full document for local use (double-click to open)
  chad_core_lab_preview.html    standalone with three.js inlined (offline screenshot testing)
"""

import json
import sys
from pathlib import Path

UI = Path(__file__).resolve().parent
ROOT = UI.parent
head = (UI / "chad_core_lab_head.html").read_text()
body = (UI / "chad_core_lab_body.html").read_text()
engine = (UI / "chad_core_lab_engine.js").read_text()
app = (UI / "chad_core_lab_app.js").read_text()
sys.path.insert(0, str(UI))
from ledger_rows import read_rows  # noqa: E402

ui_data = json.load(open(ROOT / "results" / "ui_data.json"))
rows, n_lines = read_rows()
ui_data["designs"] = rows                      # snapshot of results/design_ledger.jsonl (the Designs tab)
ui_data["designs_built"] = __import__("datetime").datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")
data = json.dumps(ui_data, separators=(",", ":"), default=float)
print(f"designs snapshot: {len(rows)} evaluations ({n_lines} ledger lines)")

page = head + body.replace("__UI_DATA__", data) + "<script>\n" + engine + "\n" + app + "\n</script>\n"
(UI / "chad_core_lab.html").write_text(page)
(UI / "chad_core_lab_standalone.html").write_text('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width,initial-scale=1">\n' + head + "</head>\n<body>\n" + page[len(head):] + "</body>\n</html>\n")
three = Path("/tmp/package/build/three.min.js")
if three.exists():
    inl = page.replace('<script src="https://cdnjs.cloudflare.com/ajax/libs/three.js/r128/three.min.js"></script>', "<script>" + three.read_text() + "</script>")
    (UI / "chad_core_lab_preview.html").write_text('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n' + head + "</head>\n<body>\n" + inl[len(head):] + "</body>\n</html>\n")
print(f"built: {len(page)/1e6:.2f} MB page, three.js inlined: {three.exists()}")
