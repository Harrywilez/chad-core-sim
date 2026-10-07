"""Serve the project folder and a live view of the design ledger for the Chad Core Lab.

    python3 ui/serve.py            # then open http://localhost:8765/ui/chad_core_lab_standalone.html

While an agent (or `python3 -m ccsim.evaluate …`) appends to results/design_ledger.jsonl, the
Lab's Designs tab polls /api/ledger and shows every new evaluation as it lands; with
"follow the newest evaluation" ticked it also rebuilds the 3-D view for each new design.

Endpoints:  /api/ledger?after=N   compact rows for ledger lines after N (+ the new count)
            /api/status            results/agent_status.json if an agent maintains one
            /api/design/<file>     a designs/*.json file
Everything else is static from the project root.  No caching, local use only.
"""

from __future__ import annotations

import json
import sys
import webbrowser
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "ui"))
from ledger_rows import read_rows  # noqa: E402

PORT = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1].isdigit() else 8765


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *a, **kw):
        super().__init__(*a, directory=str(ROOT), **kw)

    def log_message(self, fmt, *args):  # quiet
        if "/api/" not in (args[0] if args else ""):
            super().log_message(fmt, *args)

    def _json(self, obj, code=200):
        body = json.dumps(obj, default=float).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Cache-Control", "no-store")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def end_headers(self):
        self.send_header("Cache-Control", "no-store")
        super().end_headers()

    def do_GET(self):
        u = urlparse(self.path)
        if u.path == "/api/ledger":
            after = int(parse_qs(u.query).get("after", ["0"])[0])
            rows, n = read_rows(after)
            return self._json({"rows": rows, "count": n})
        if u.path == "/api/status":
            f = ROOT / "results" / "agent_status.json"
            return self._json(json.load(open(f)) if f.exists() else {})
        if u.path.startswith("/api/design/"):
            f = ROOT / "designs" / Path(u.path.split("/api/design/", 1)[1]).name
            return self._json(json.load(open(f))) if f.exists() else self._json({"error": "no such design"}, 404)
        return super().do_GET()


if __name__ == "__main__":
    url = f"http://localhost:{PORT}/ui/chad_core_lab_standalone.html"
    print(f"Chad Core Lab (live ledger): {url}\nproject root: {ROOT}\nCtrl-C to stop")
    if "--no-browser" not in sys.argv:
        try:
            webbrowser.open(url)
        except Exception:
            pass
    ThreadingHTTPServer(("127.0.0.1", PORT), Handler).serve_forever()
