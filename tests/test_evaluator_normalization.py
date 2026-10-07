"""Evaluator v2.3 normalisation tests (fast settings; numba or numpy).

    python3 tests/test_evaluator_normalization.py

1. A uniform rescaling of every relative current is physically a no-op and must leave score,
   tau_c and buildability unchanged (events 122/126/139 showed a 1/k change before v2.3).
2. Reordering components must not change the current normalisation or buildability.
3. A zero-current component is disabled: same result as the design without it.
4. Joule power must reflect relative amplitudes (a 0.3-relative coil dissipates 0.09× per unit length).
"""

from __future__ import annotations

import copy
import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ccsim.design import build  # noqa: E402
from ccsim.evaluate import evaluate  # noqa: E402

FAST = dict(transits=4, n_particles=16, grid_n=9, max_transits=4, verbose=False, log=False)
BASE = {"name": "norm-test", "wire_radius": 0.01,
        "components": [{"family": "solenoid", "params": {"radius": 0.45, "length": 1.1, "turns": 8}, "current": 1.0},
                       {"family": "circle", "params": {"radius": 0.7, "axis": "z"}, "translate": [0, 0, 0.7], "current": 0.3}]}


def close(a, b, rel=1e-6):
    return abs(a - b) <= rel * max(abs(a), abs(b), 1e-12)


def test_uniform_scaling_invariance():
    r1 = evaluate(BASE, **FAST)
    d2 = copy.deepcopy(BASE)
    for c in d2["components"]:
        c["current"] *= 1.56108445
    r2 = evaluate(d2, **FAST)
    # the normalisation is exactly invariant …
    assert close(r1["buildability"]["factor"], r2["buildability"]["factor"], 1e-9)
    assert close(r1["buildability"]["I_max_A"], r2["buildability"]["I_max_A"], 1e-9), (r1["buildability"]["I_max_A"], r2["buildability"]["I_max_A"])
    assert close(r1["field"]["B_rms_per_364kA_T"], r2["field"]["B_rms_per_364kA_T"], 1e-9)
    assert close(r1["field"]["current_A"], r2["field"]["current_A"] * 1.56108445, 1e-9)       # per-unit-relative current scales as 1/k
    assert close(r1["field"]["current_first_winding_A"], r2["field"]["current_first_winding_A"], 1e-9)  # same physical current
    # … and the orbits are the same field to rounding: all but at most one chaotic long-lived marker agree
    l1 = [x or float("inf") for x in r1["particles"]["loss_times_transits"]]
    l2 = [x or float("inf") for x in r2["particles"]["loss_times_transits"]]
    differing = sum(1 for a, b in zip(l1, l2) if not (a == b or close(a, b, 1e-3)))
    assert differing <= 1, (differing, l1, l2)


def test_reorder_invariance():
    r1 = evaluate(BASE, **FAST)
    d2 = copy.deepcopy(BASE)
    d2["components"] = d2["components"][::-1]
    r2 = evaluate(d2, **FAST)
    assert close(r1["buildability"]["I_max_A"], r2["buildability"]["I_max_A"], 1e-6), (r1["buildability"]["I_max_A"], r2["buildability"]["I_max_A"])
    assert close(r1["buildability"]["factor"], r2["buildability"]["factor"], 1e-6)
    assert close(r1["field"]["B_rms_roi_T"], r2["field"]["B_rms_roi_T"], 1e-9)
    # the physical current of the unit-relative component is the same whichever component comes first
    assert close(r1["field"]["current_A"], r2["field"]["current_A"], 1e-6)


def test_zero_current_component_is_disabled():
    d0 = {"name": "solo", "wire_radius": 0.01, "components": [BASE["components"][0]]}
    dz = copy.deepcopy(BASE)
    dz["components"][1]["current"] = 0.0
    assert len(build(dz).windings) == len(build(d0).windings)
    r0, rz = evaluate(d0, **FAST), evaluate(dz, **FAST)
    assert close(r0["score"], rz["score"], 1e-9) and close(r0["buildability"]["I_max_A"], rz["buildability"]["I_max_A"], 1e-9)


def test_joule_power_uses_relative_amplitudes():
    w = build(BASE).with_current(1000.0)
    I = [abs(x.current_A) for x in w.windings]
    assert close(I[0], 1000.0) and close(I[1], 300.0), I


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("ok", name)
