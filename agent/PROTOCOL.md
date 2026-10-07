# Design-agent protocol

You are optimising magnetic-confinement coil assemblies in `ccsim`. Everything you need is
one command away and every number comes from a simulation of the design you propose.

## The objective

`score = τ_c × buildability`,  `τ_c = √(τ_all · τ_pass)`

* **τ_all** — confinement time in wall transits of 15 keV D⁺ markers seeded uniformly in the
  vessel (r ∈ [0.18, 0.58] of a unit ball, ≥ 0.03 from any conductor), isotropic pitch,
  full-orbit Boris on the exact vacuum field at **B_rms = 2 T** in the region of interest
  (the current is set to make that so — every design is judged at the same field strength).
  τ = ∫S dt plus a fitted exponential tail capped at one more run-length; 40 transits,
  extended to 80 and 120 while ≥ 10 % of the markers are still confined and the tail is a
  large share of τ (so plateaus are resolved rather than extrapolated).
* **τ_pass** — the same for the *passing* markers: those inside the loss cone of their own
  field line (ξ0² ≥ 1 − 1/R_line) plus every marker on a closed line. Mirror trapping cannot
  hold these; only closed field lines / flux surfaces can. A plasma needs them held (they are
  refilled by collisions), so a trap that keeps only trapped pitch angles — a mirror, a
  dipole-like fat torus — is not rewarded as if it held everything. That is why the fat
  single knot that scored 300 under the first version of the score scores 6 now.
* **buildability** = min(1, J_REBCO / J), J = I_max/(π a²) with a = min(design wire radius,
  0.45 × clearance) — the same conductor radius the field model uses — and I_max including
  the largest relative component current; J_REBCO = 10⁹ A/m². Declaring a thinner wire than
  the clearance allows now only hurts. A design whose conductors touch (clearance < 2 × wire
  radius) or leave the unit ball is invalid (score 0).
* **score_plasma** (only with `--surfaces`, toroidal designs, ~2 min): score × a rung-3
  factor — a toroidal-field configuration without enough rotational transform to short the
  polarization is capped at τ_pol; one with ≥ 20 % nested surfaces and
  |ι_edge| ≥ max(0.02, (v_d,i/v_th,e)·πR/a) — the electrons carry the shorting current, so
  a few hundredths of transform suffice — is exempt. Since v2.2 the cap applies both to
  closed-line designs and to *leaky* toroidal solenoids (mirror ratio < 1.5: nothing is
  mirror-trapped, the markers are only held by a slow drift of the lines to the wall);
  mirror designs (R ≥ 1.5) get no factor because their passing-particle loss is already in
  τ_pass. Use it to rank finalists.

Secondary numbers you get for free and should reason with: S4/S20/S120 for all and for
passing markers, the closed-field-line fraction, mirror ratio and adiabatic prediction, the
current and Joule power, and — with `--surfaces` — flux-surface fraction, last closed
surface radius, ι and magnetic well.

The score's first version (τ_all with an uncapped tail, a wire-radius term that could be
gamed) is what the first agent run optimised; its report is in agent/REPORT_run1.md and
its leaderboard was re-scored under this version in agent/leaderboard.md.

## The commands

```bash
python3 -m ccsim.evaluate designs/<name>.json                 # ~10–30 s, prints one line + JSON
python3 -m ccsim.evaluate designs/<name>.json --surfaces      # + Poincaré / ι (toroidal only, +1–3 min)
python3 -c "from ccsim.design import FAMILIES; import json; print(json.dumps(FAMILIES, indent=1))"
python3 agent/optimize.py designs/<seed>.json --iters 30      # baseline random-mutation hill climb
tail -n 20 results/design_ledger.jsonl                        # every evaluation ever made (append-only)
```

## Geometry you must know about (session 8)

The meshed-torus families (`hopfmesh`, `hopfmirror` with 2 revolutions) use a *zone weave*: the
two strands interlock only where they are far enough apart in the surface for an over/under
swap; where the modulated turns run close the over/under is frozen. Consequence: **δη ≲ 0.015
weaves, δη ≳ 0.03 nests** (strand A on the outer torus, B on the inner one, no interlock, a
coaxial gap field, no closed lines for the chiral design). The clearance reported by
`validate()` is exact; the ~0.02 same-strand turn spacing of a 24-circuit torus is the ceiling.
Designs evaluated before 2026-09-05 used a geometry whose strands touched — treat those ledger
rows as history, not as data.

## The design format

```json
{"name": "…", "wire_radius": 0.005,
 "components": [
   {"family": "precess", "params": {"inward": 0.22, "rotation_deg": 24, "slip_deg": 3},
    "scale": 0.45, "reflect": null, "rotate": {"axis": [1,0,0], "deg": 90},
    "translate": [0.5, 0, 0], "current": 1.0}
 ]}
```
Transform order: scale → reflect (x|y|z) → rotate → translate. `current` is relative (sign =
direction); multi-circuit families (hopfmirror, hopfmesh, yinyang, picket, link) carry their
own internal signs (`sense`). The `cone` family is a conical helix ("tornado coil": base at
−height/2, tip at +height/2; reflect z to point the tip down) for plug / nozzle experiments. All conductors must stay inside r ≤ 1; the vessel wall is at
r = 0.82 and particles are seeded at r ≤ 0.58, so a component's active region should sit
inside ~0.8. Point budget: ≤ 60 000 conductor points per design.

## How to work

1. Start from the seeds in `designs/` and their ledger lines. Read PLAYBOOK.md §5 and RUNGS.md
   §0 first — they tell you which physics moves the number (closed lines beat mirror ratio;
   the drift tax scales with ρ_L/L_B; fatter tori win; an achiral pair has no transform).
2. Propose designs from a **stated hypothesis** ("moving the two cores apart lengthens L_B,
   so S20 should rise"), evaluate, and write one line in `agent/notes.md`: hypothesis →
   result → what it changed about your model. Falsified hypotheses are as valuable as
   confirmed ones; record them.
3. Keep a leaderboard (`agent/leaderboard.md`): score, τ_c, buildability, design file.
4. Prefer changes that a physicist would make (clearance, core separation, current sense,
   handedness, layering a chiral weave over a mirror) over blind parameter noise; use
   `agent/optimize.py` for local polishing once a promising basin is found.
5. Respect the budget you were given (number of evaluations or wall-clock). Each
   evaluation is ~10–30 s; a `--surfaces` evaluation is ~2 min.
6. Never edit `ccsim/` to make a number better. If you find a bug or an unfair term in the
   score, write it down in `agent/notes.md` and keep going with the score as it is.

## What the score cannot see (say so in your report)

Collisions, plasma pressure and the self-consistent electric field (RUNGS.md §3: a
toroidal field with no transform loses its plasma by polarization even when single
particles look fine — use `--surfaces` and ι to check toroidal designs), heating and
fuelling, wall interaction, neutronics. The score rewards single-particle vacuum-field
confinement per unit of conductor stress; treat it as the first filter, not the verdict.
