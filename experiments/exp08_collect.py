"""Merge exp08 parts (and the exp07 rows) into results/exp08_playbook.json and the UI ledger."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
R = ROOT / "results"

merged = {}
for p in sorted(R.glob("exp08_*.json")):
    if p.name in ("exp08_playbook.json",):
        continue
    merged.update(json.load(open(p)))


def ui_spec(name: str, meta: dict) -> dict:
    fam = meta.get("family")
    if fam == "hopf-mirror":
        return {"family": "hopfmirror", "params": {"eta": meta["eta"], "circuits": meta["circuits"], "revs": meta.get("revolutions", 1), "mode": meta["mode"],
                                                   "deta": meta.get("delta_eta", 0.025), "sense": "same" if meta.get("sense", 1) > 0 else "opposed"}, "layout": "single"}
    if fam == "hopf-torus":
        return {"family": "hopftorus", "params": {"eta": meta["eta"], "circuits": meta["circuits"], "revs": meta["revolutions"]}, "layout": "single"}
    if fam == "hopf-continued":
        return {"family": "hopfcont", "params": {"eta0": 0.82, "eta1": meta["eta1"], "sweep": meta["sweep"], "circuits": meta["circuits"]}, "layout": "single"}
    if fam == "hopf-drift":
        return {"family": "hopfdrift", "params": {"variant": meta["variant"], "circuits": meta["circuits"]}, "layout": "single"}
    if fam == "baseball":
        return {"family": "baseball", "params": {"amp": meta["amplitude_deg"], "turns": meta["turns"]}, "layout": "single"}
    if fam == "yin-yang":
        return {"family": "yinyang", "params": {"turns": meta["turns"]}, "layout": "single"}
    if fam == "picket":
        return {"family": "picket", "params": {"rings": meta["rings"]}, "layout": "single"}
    if fam == "sphere":
        return {"family": "sphere", "params": {"turns": meta["turns"]}, "layout": "single"}
    if fam == "link":
        return {"family": "link", "params": {"sense": "opposed" if meta.get("sense", 1) < 0 else "same"}, "layout": "single"}
    if fam == "assembly":
        core = meta["core"]
        if core == "precess":
            return {"family": "precess", "params": {"inward": meta.get("inward", 0.22), "rotation": meta.get("rotation_deg", 24), "slip": 3},
                    "layout": meta["layout"], "corescale": meta.get("core_scale", 0.45)}
        if core == "baseball":
            return {"family": "baseball", "params": {}, "layout": meta["layout"], "corescale": 0.45}
        if core == "yin-yang":
            return {"family": "yinyang", "params": {"sense": "opposed" if meta.get("inner_sign", 1) < 0 else "same"}, "layout": meta["layout"], "corescale": meta.get("core_scale", 0.5)}
    return {}


rows = []
for name, rec in merged.items():
    t, f, r = rec["topology"], rec["field"], rec["reactor"]
    rows.append({"name": name, "closed": t["frac_closed"], "wall": t["frac_wall_both"], "Rmir": t["mirror_ratio_p90"], "pred": t["predicted_adiabatic_retention"],
                 "S4": r["survival"]["S4"], "S8": r["survival"]["S8"], "S20": r["survival"]["S20"], "median_loss_transits": r["median_loss_transits"],
                 "core": r["core_fraction"], "core0": r["core_fraction_initial"], "brms": f["B_rms_roi_T"], "bmax_rms": f["B_max_over_rms"],
                 "clearance": rec["clearance"], "active_length": rec["active_length"], "joule_MW": r["joule_power_W"] / 1e6, "energy_drift": r["energy_drift_rel_max"],
                 "ui": ui_spec(name, rec["meta"]), "meta": rec["meta"]})
rows.sort(key=lambda x: -x["S20"])
json.dump(rows, open(R / "exp08_playbook.json", "w"), indent=1, default=float)

ui = json.load(open(R / "ui_data.json"))
ui["playbook"] = [{k: v for k, v in row.items() if k != "meta"} for row in rows]
# rung 1: flux-surface scan (exp10) → ledger rows
import glob
surf = {}
for pth in sorted(glob.glob(str(R / "exp10_*.json"))):
    surf.update(json.load(open(pth)))
srows = []
for name, r in surf.items():
    m = r.get("meta", {})
    fam = m.get("family", "")
    if fam == "hopf-mirror-chiral":
        uispec = {"family": "hopfmesh", "params": {"eta": m["eta"], "circuits": m["circuits"], "eps": m["eps"], "periods": m["periods"], "deta": m.get("delta_eta", 0.03), "hmode": m.get("helical_mode", "toroidal"), "sense": "opposed"}, "layout": "single"}
    elif fam == "hopf-mirror-twisted" and m.get("ellipticity", 0) == 0 and m.get("helical_axis", 0) == 0:
        uispec = {"family": "hopfmesh", "params": {"eta": m["eta"], "circuits": m["circuits"], "eps": 0, "periods": m["periods"], "deta": 0.03, "hmode": "toroidal", "sense": "opposed"}, "layout": "single"}
    else:
        uispec = {}
    srows.append({"name": name, "family": fam, "surf": r["frac_surface"], "island": r["frac_island_chaotic"], "open": r["frac_open"], "r_out": r["outermost_surface_r"],
                  "iota_axis": r["iota_axis"], "iota_edge": r["iota_edge"], "well": r["well_depth"], "axis_R": r["axis_R"], "clearance": r["clearance"], "ui": uispec})
srows.sort(key=lambda x: -(abs(x["iota_edge"]) if x["iota_edge"] == x["iota_edge"] else 0))
ui["surfaces"] = srows
tr = {}
for pth in sorted(glob.glob(str(R / "exp11_*.json"))):
    tr.update(json.load(open(pth)))
ui["transport"] = [{"name": k, **v["survival"], "nu": v["nu_per_transit"], "r_seed": v["r_seed"]} for k, v in tr.items()]
v1 = json.load(open(R / "exp01_verify.json"))
if "fields_sphere_winding" in v1:
    ui["verify"]["sphere"] = v1["fields_sphere_winding"]
json.dump(ui, open(R / "ui_data.json", "w"), default=float)
print(f"{len(rows)} configurations")
if "--print" in sys.argv:
    print(f"{'configuration':55s} {'closed':>6s} {'Rmir':>5s} {'pred':>5s} {'S4':>5s} {'S8':>5s} {'S20':>5s} {'core':>5s} {'Brms':>6s} {'Bmax/rms':>8s} {'MW':>6s}")
    for x in rows:
        print(f"{x['name']:55s} {x['closed']:6.2f} {x['Rmir']:5.2f} {x['pred']:5.2f} {x['S4']:5.2f} {x['S8']:5.2f} {x['S20']:5.2f} {x['core']:5.2f} {x['brms']:6.2f} {x['bmax_rms']:8.1f} {x['joule_MW']:6.1f}")
