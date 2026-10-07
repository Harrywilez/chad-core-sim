"""README.md → a self-contained HTML page (figures embedded, MathJax for the equations).

    python3 experiments/export_hopf_family_html.py            # writes exports/hopf_family/hopf_family.html
                                                              # and exports/hopf_family/hopf_family_artifact.html
The artifact variant has no <html>/<head>/<body> skeleton (the Artifact tool wraps it).
"""

from __future__ import annotations

import base64
import re
import sys
from pathlib import Path

import markdown

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "exports" / "hopf_family"

CSS = """
:root{--bg:#f6f7f9;--ink:#1b2430;--muted:#5d6b7a;--rule:#d9dee5;--accent:#2a73d9;--counter:#d94d33;--code:#eceff3;--figbg:#ffffff}
@media (prefers-color-scheme: dark){:root:not([data-theme="light"]){--bg:#121820;--ink:#e4e9ef;--muted:#9fabb8;--rule:#283442;--accent:#7fb0f5;--counter:#f28a70;--code:#1b2530;--figbg:#f6f7f9}}
:root[data-theme="dark"]{--bg:#121820;--ink:#e4e9ef;--muted:#9fabb8;--rule:#283442;--accent:#7fb0f5;--counter:#f28a70;--code:#1b2530;--figbg:#f6f7f9}
html{background:var(--bg)}
body{background:var(--bg);color:var(--ink);margin:0;font-family:"IBM Plex Sans",-apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;font-size:17px;line-height:1.6;-webkit-font-smoothing:antialiased}
main{max-width:1000px;margin:0 auto;padding:40px 24px 96px}
.lede{max-width:none}
.eyebrow{font-family:"IBM Plex Mono",ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.78em;letter-spacing:.08em;text-transform:uppercase;color:var(--muted);margin:0 0 14px}
.eyebrow .a{color:var(--accent)} .eyebrow .b{color:var(--counter)}
h1{font-family:"IBM Plex Serif",Georgia,"Times New Roman",serif;font-weight:600;font-size:2.05em;line-height:1.15;letter-spacing:-.01em;margin:0 0 .5em;text-wrap:balance;max-width:26ch}
h2{font-family:"IBM Plex Serif",Georgia,serif;font-weight:600;font-size:1.45em;line-height:1.25;margin:2.2em 0 .6em;padding-top:.9em;border-top:1px solid var(--rule);text-wrap:balance;max-width:72ch}
h3{font-family:"IBM Plex Serif",Georgia,serif;font-weight:600;font-size:1.12em;margin:1.6em 0 .4em;text-wrap:balance;max-width:72ch}
p,ul,ol{max-width:72ch}
p{margin:0 0 1em}
ul,ol{margin:0 0 1em;padding-left:1.4em}
li{margin:.25em 0}
a{color:var(--accent);text-decoration-thickness:1px;text-underline-offset:2px}
strong{font-weight:600}
code{font-family:"IBM Plex Mono",ui-monospace,SFMono-Regular,Menlo,monospace;font-size:.86em;background:var(--code);padding:.08em .35em;border-radius:3px}
pre{max-width:72ch;margin:1em 0;background:var(--code);padding:14px 16px;border-radius:4px;overflow-x:auto;font-size:.84em;line-height:1.5}
pre code{background:none;padding:0;font-size:1em}
table{border-collapse:collapse;margin:1.2em 0;font-size:.88em;display:block;overflow-x:auto;max-width:100%;width:max-content;font-variant-numeric:tabular-nums}
th,td{border-bottom:1px solid var(--rule);padding:6px 12px;text-align:left;vertical-align:top}
th{font-family:"IBM Plex Mono",ui-monospace,Menlo,monospace;font-weight:500;font-size:.82em;letter-spacing:.04em;text-transform:uppercase;color:var(--muted);border-bottom:2px solid var(--rule)}
td:first-child{color:var(--muted)}
img{max-width:100%;height:auto;display:block;margin:1.6em 0;background:var(--figbg);border:1px solid var(--rule);border-radius:3px}
mjx-container[display="true"]{overflow-x:auto;overflow-y:hidden;max-width:72ch;margin:.6em 0 !important}
.meta{font-family:"IBM Plex Mono",ui-monospace,Menlo,monospace;font-size:.78em;color:var(--muted);max-width:72ch;margin:2.5em 0 0;padding-top:1em;border-top:1px solid var(--rule)}
@media (max-width:640px){body{font-size:16px} main{padding:24px 16px 64px} h1{font-size:1.7em}}
@media (prefers-reduced-motion: reduce){*{scroll-behavior:auto}}
"""

MATHJAX = """
<script>
window.MathJax = {tex: {inlineMath: [['$','$'], ['\\\\(','\\\\)']], displayMath: [['$$','$$']], processEscapes: true},
                  svg: {fontCache: 'global'}, options: {skipHtmlTags: ['script','noscript','style','textarea','pre','code']}};
</script>
<script async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-svg.js"></script>
"""


def protect_math(text: str):
    """Replace $$…$$ and $…$ spans by placeholders so the markdown converter leaves them alone."""
    store = []

    def keep(m):
        store.append(m.group(0))
        return f"@@MATH{len(store) - 1}@@"

    text = re.sub(r"\$\$.*?\$\$", keep, text, flags=re.S)
    text = re.sub(r"(?<!\\)\$(?!\s)(.+?)(?<!\s)\$", keep, text)
    return text, store


def restore_math(html: str, store: list) -> str:
    def back(m):
        s = store[int(m.group(1))]
        return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return re.sub("@@MATH(\\d+)@@", back, html)


def embed_images(html: str) -> str:
    def rep(m):
        src = m.group(1)
        p = OUT / src
        if p.exists():
            data = base64.b64encode(p.read_bytes()).decode("ascii")
            return f'src="data:image/png;base64,{data}"'
        return m.group(0)
    return re.sub(r'src="(figures/[^"]+)"', rep, html)


def main() -> None:
    md = (OUT / "README.md").read_text(encoding="utf-8")
    text, store = protect_math(md)
    body = markdown.markdown(text, extensions=["tables", "fenced_code"])
    body = restore_math(body, store)
    body = embed_images(body)
    title = "The Hopf Family"
    desc = "The Hopf-drift winding and the two versions of the Hopf torus + mirror image: equations, files, fields."
    fonts = '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Sans:ital,wght@0,400;0,600;1,400&family=IBM+Plex+Serif:wght@600&family=IBM+Plex+Mono:wght@400;500&display=swap">'
    head = f"<title>{title}</title>\n{fonts}\n<style>{CSS}</style>\n{MATHJAX}\n"
    body = body.replace("<h1>", '<p class="eyebrow">Chad Core Sim · exports/hopf_family · <span class="a">strand A blue</span> · <span class="b">strand B red</span></p>\n<h1>', 1)
    inner = f'<main>\n<div class="lede">{body}</div>\n<p class="meta">Generated from exports/hopf_family/README.md by experiments/export_hopf_family_html.py; every number is in MANIFEST.json.</p>\n</main>\n'
    (OUT / "hopf_family_artifact.html").write_text(head + inner, encoding="utf-8")
    # the local file inlines MathJax (2.1 MB) so it renders offline; the artifact loads it from cdnjs (CSP allow-list)
    mj = OUT / "figures" / ".mathjax-tex-svg-3.2.2.js"
    head_local = head.replace('<script async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/3.2.2/es5/tex-svg.js"></script>',
                              "<script>" + mj.read_text(encoding="utf-8") + "</script>") if mj.exists() else head
    full = ('<!doctype html>\n<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">\n'
            + head_local + "</head><body>\n" + inner + "</body></html>\n")
    (OUT / "hopf_family.html").write_text(full, encoding="utf-8")
    print("wrote", OUT / "hopf_family.html", len(full) // 1024, "kB")


if __name__ == "__main__":
    main()
