"""run-2 helper: write a design JSON with wire_radius = 0.45 x clearance (the value that
maximises a = min(wire_radius, 0.45*clearance) under score v2 while staying valid), and/or
report geometry for a spec without writing.

Usage:
  python3 agent/mk2.py write designs/out.json '<json spec>'
  python3 agent/mk2.py scan '<json list of specs>'
"""
import json
import sys

sys.path.insert(0, '/home/claude/chadcore-sim')
from ccsim.design import validate  # noqa: E402


def fit(spec):
    spec = dict(spec)
    spec.setdefault("wire_radius", 0.0005)
    v = validate(spec)
    spec["wire_radius"] = round(0.45 * v["clearance"], 6)
    return spec, validate(spec)


def main():
    mode = sys.argv[1]
    if mode == "write":
        out, spec = sys.argv[2], json.loads(sys.argv[3])
        spec, v = fit(spec)
        json.dump(spec, open(out, "w"), indent=1)
        print(f"{out}: clr={v['clearance']:.4f} wr={spec['wire_radius']:.5f} rmax={v['r_max']:.3f} "
              f"pts={v['conductor_points']} valid={v['valid']} {v['problems']}")
    else:
        specs = json.loads(sys.argv[2])
        for spec in specs:
            s, v = fit(spec)
            print(f"{spec.get('name','?')[:46]:46s} clr={v['clearance']:.4f} wr={s['wire_radius']:.5f} "
                  f"rmax={v['r_max']:.3f} pts={v['conductor_points']:6d} valid={v['valid']} {v['problems']}")


main()
