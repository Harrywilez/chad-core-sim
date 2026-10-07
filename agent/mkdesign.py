"""Write a design JSON with wire_radius set to clearance/3 (the value that makes the
score's clearance term exactly 1 — see notes.md F1).  Usage:
  python3 agent/mkdesign.py designs/out.json '<json with components/name>'
"""
import json, sys
sys.path.insert(0, '/home/claude/chadcore-sim')
from ccsim.design import validate

out, spec = sys.argv[1], json.loads(sys.argv[2])
spec.setdefault("wire_radius", 0.001)
v = validate(spec)
spec["wire_radius"] = round(v["clearance"] / 3.0, 5)
v2 = validate(spec)
json.dump(spec, open(out, "w"), indent=1)
print(f"{out}: clr={v['clearance']:.4f} wr={spec['wire_radius']:.5f} rmax={v['r_max']:.3f} pts={v['conductor_points']} valid={v2['valid']} {v2['problems']}")
