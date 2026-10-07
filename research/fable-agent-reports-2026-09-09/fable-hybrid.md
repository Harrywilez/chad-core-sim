# fable-hybrid — phase-2 report (2026-09-09, 05:03–06:05 UTC)

Direction: (A) plug the last losses of the two leading families with physically motivated additions; (B) a placement-aware
structural search with `agent/composition.propose()` seeded from the leaders. 47 screening evaluations + 2 confirmations at
480 markers (49 of the 60 + 3 budget, 62 min wall-clock); notes with every batch in `notes/fable-hybrid.md`; generators `designs/fable-hybrid/_gen_p1.py`, `_gen_p2.py`,
search driver `designs/fable-hybrid/_search.py` (round logs `_search_round{1,2,3}.json`, `_search_r*_*.log`).

**Headline: nothing I added beats its host beyond the ±5–10 % noise.** The magnetosphere's polar null cannot be plugged by
polar coils (they relabel wall losses as coil losses, net zero — confirmed at 480 markers with identical retained counts), a
toroidal field on the magnetosphere costs > 80 % of τ at any current, and the structural search (144 proposals, 3 rounds,
two parents) produced 0 improvements: every operator in `composition.py` acts on one component and breaks the up/down
symmetry the leaders depend on; heterogeneous adds cost 10–90 %.

## (a) Top 5 of my evaluations (fixed 120 T, 160 markers unless noted; host = the design the addition was made to)

| # | design (file in `designs/fable-hybrid/`) | score | τ_all / τ_pass | build | closed | S120 | host score |
|---|---|---|---|---|---|---|---|
| 1 | `b1_32_parameter.json` — fable-rings' belt leader with the top belt ring at R 0.923 (a 2.5 % tweak) | 210.6 | 210 / 211 | 1.00 | 1.00 | 0.86 | 210.3 (`designs/fable-rings/c1_belt_R45n2_p90z25_i60.json`) |
| 2 | `b1_24_parameter.json` — same, bottom belt ring at R 0.877 | 209.9 | 209 / 211 | 1.00 | 1.00 | 0.86 | 210.3 |
| 3 | `b2_01_current.json` — b1.32 with the top belt at −0.586 | 208.2 | 207 / 210 | 1.00 | 1.00 | 0.84 | 210.6 |
| 4 | `p1_cone_m050.json` — magnetosphere R0.40 n2 i50 + polar cones z ±0.75 tip-in, I −0.05 | 180.5 (185.7 @480 m) | 206 / 212 (212 / 219 @480) | 0.86 | 0.97 | 0.83 (0.86) | 177.5 (180.5 @480) `p1_host.json` = `designs/fable-rings/b3_vfs_R40_n2_i50.json` |
| 5 | `p1_cone_m010.json` — same host, cones I −0.01 | 179.3 | 207 / 209 | 0.86 | 1.00 | 0.84 | 177.5 |

Rows 1–3 are the leader ± noise (single-component parameter/current moves of the search); rows 4–5 are the polar plug
(neutral, see below). Nothing here is a new leader; the room's leader remains fable-rings' belt design.

## (b) Plug findings (Part A, 25 evaluations)

1. **Polar cones on the magnetosphere — H falsified both signs.** Calibrated first with Biot–Savart on the axis: the host's
   on-axis null sits at z = 0.677 with a gradient ≈ 1 T/m (non-adiabatic radius ≈ 0.15); a cone (base .30, tip .06, h .30,
   6 turns) at z ±0.75 tip-in gives 6.5–7.6 × the host's centre field per unit current at the null, so the task's 0.05–0.3 range
   dominates the whole polar region and the useful range is 0.002–0.01. Results (wall / wire / retained of 160, host
   10/16/134): I = −0.002…−0.01 → 3–5 / 20–23 / 133–136; +0.002…+0.01 → 5–7 / 19–21 / 134; ±0.05 → 4/23/133 (polar mirror,
   R_mirror 8) and 13/29/118 (separatrix opened, closed 0.92). Retained counts are flat; the antiparallel cone merely
   intercepts the escaping polar column at its tip (wall → "wire"). Confirmed at 480 markers: cone −0.05 retained 413
   (wall 12, wire 55) vs host 412 (wall 28, wire 40); its +3 % score is τ_pass bookkeeping (23 markers become mirror-trapped
   and leave the passing set). Physics: the null is topologically protected on the axis (B there is purely axial and must
   change sign between the ring and the opposing field); an axisymmetric coil can only move it, and the non-adiabatic zone
   moves with it. Polar rings (R .40, z ±.86, ±0.02) behave the same (−0.02 flat, +0.02 −4 %).
2. **Toroidal field on the magnetosphere — the cost curve is a cliff.** A nested `hopfmirror` η 1.0 δη 0.2 (N 4, clearance
   0.040 so the global wire radius barely changes) at 5 / 10 / 20 / 40 % of the ring-turn current: τ 207 → 38 / 32 / 24 / 6
   (score 27 / 22 / 16.5 / 4; closed 0.85 / 0.62 / 0.46 / 0.08; S20 0.28 at 5 %). Unaliased ι (`tools/iota.py`, axis = the ring)
   at 5 %: −219 near the ring, −35 at 0.22, −4.6 at 0.36 — transform is never the constraint, confinement is. Ripple control:
   N 24 woven at 10 % keeps τ 113 (closed 0.73) — so ~75 % of the N 4 loss is the 8-strand ripple destroying the drift
   surfaces (p_φ is no longer conserved), the rest is strand interception in the hole (ρ 0.03–0.06) and at the tube top
   (z ≈ 0.45); but N 24 has clearance 0.016 → wire at the cap → build 0.21 → score 24. Either way the add-on is dead on this
   evaluator; the argument for the dipole plasma not polarizing (brief §7) has to carry the plasma case.
3. **Levitron + external PF / cusp rings — no plug.** Opposing rings (0.70, ±0.60) at −0.3 / −0.5 / −0.8: 91.0 / 91.7 / 88.3 vs
   host 88.5 (same-batch control) — noise; the closed fraction falls (0.44 → 0.31) while τ_all rises 84 → 88. Cusp pair ±0.5:
   69.6. The Levitron's losses are the strand curtain (every line that crosses the η 0.9 tube surface dies in a few bounces;
   the polar column is lost promptly) and ripple transport of the lines inside the tube; nothing outside the wall reaches them.
4. **Split / wide ring stacks hurt both hosts.** Magnetosphere rings at z ±0.08 instead of ±0.02: 142 (wire losses 16 → 43).
   Levitron ring 2.0 → 2 × 1.0 at z ±0.06 / ±0.10: 53 / 80 vs 88. A wider stack is a bigger wire target and a weaker field slot.

## (c) Structural search (Part B, 23 evaluations)

Driver: `propose()` with `agent/optimize.mutate` as the parameter mutator, sigma 0.15 (0.30 in round 3), structural rate 0.65,
a 13-donor library of sensibly placed components (polar / belt / hole / cap rings, extra stack turn, polar cone, `picket` PF
pair, cusp pair and shell, `sphere`, `torushelix`, nested N 4 toroidal solenoid, the η 0.9 Levitron torus); every proposal is
repaired (wire = 0.45 × clearance), validated, and pre-screened on a 15³ grid (closed fraction of 48 shell seeds + buildability
estimate, ~8 s) — "plausible" = closed ≥ 0.5 × parent and build ≥ 0.3. Rounds 2–3 add a symmetrize repair (mirror any added
off-midplane component) and reserve slots for heterogeneous candidates.

| round | parent | proposals | valid | plausible | submitted | improved | best child |
|---|---|---|---|---|---|---|---|
| 1 | belt leader 210.3 | 40 | 29 (6 outside the ball, 5 duplicates) | 21 (3 heterogeneous) | 8 | 0 (2 within noise) | 210.6 (belt ring R .923) |
| 1 (Levitron) | L_e0.90_N24_r0.400 88.5 | 24 | 14 (everything placed outside collides with the tube) | 7 (1 heterogeneous) | 2 | 0 | 86.1 (single off-axis cone add) |
| 2 | b1.32 210.6 | 40 | 30 | 22 (5 heterogeneous) | 8 | 0 (1 within noise) | 208.2 (belt −0.586); sphere add 189, ring-pair adds 94–102 |
| 3 | b2.01 208.2 (sigma 0.30) | 40 | 29 | 17 (4 heterogeneous) | 4 | 0 | 193.9 (third belt ring at the equator, closed 0.86); picket pair add 156.5 |

Overall: 144 proposals → 102 valid (71 %) → 67 plausible (47 %) → 22 submitted → 0 improved. Compared with the phase-1
optimizer run (7 invalid of 12, no heterogeneous valid candidate), the placed donor library + repair + pre-screen lifts the
valid rate to 70 % and does produce heterogeneous valid candidates (sphere, torushelix, picket, ring pairs) — but every one
of them scores below its parent (sphere −10 %, ring pairs −50 %, torushelix −90 %), and the homogeneous moves that survive
are single-ring tweaks that break the up/down symmetry (−1 to −25 %). What the search found is that the leader sits on a
symmetric ridge: any one-sided change costs, and the operator set has no symmetric-pair move except my mirror repair.

## (d) Lessons

1. The magnetosphere's remaining losses are (i) seeds born within ~0.05–0.13 of the coil (wire, ~7–8 %) and (ii) the polar
   column through the on-axis null (~6–14 %); (ii) can be reduced by reshaping the separatrix (fable-rings' belt) but not by
   plugging the null, which is topologically protected on the axis for any axisymmetric coil set.
2. A converging polar coil ("tornado") is a mirror plug only above ~0.05 × ring current, and at that strength it opens the
   separatrix (+ sign) or intercepts the column on its own tip (− sign). Retained counts, not wall/wire labels, are the metric.
3. Any toroidal-field winding that fits the magnetosphere either ripples (N ≤ 6: τ −80 %) or thins the global wire
   (N 24: build 0.2); the closed poloidal lines are also cut by the tube surface. Toroidal field and levitated dipole don't mix
   in this geometry.
4. Two rings are only better than one at the chord-pitch spacing: any wider slot (±0.06–0.10) multiplies wire hits.
5. The Levitron is a dead end for additions: its torus fills the ball (outer edge r = 1.0), so nothing fits outside it, and its
   losses are inside the tube (strand curtain, ripple).
6. Structural search on a symmetric leader needs symmetric operators; `composition.propose()` has none, and its jitter pushes
   any donor at ρ ≳ 0.9 out of the ball. A coarse-grid pre-screen (closed × build) ranks candidates well enough to avoid wasting
   evaluations on the junk (its one miss: torushelix closed 0.69 on the 15³ grid vs 0.45 in the evaluator).
7. Cross-check the noise: the same host evaluated three times (vm 178.5, vm 177.5, 480 m 180.5) — differences < 2 %.

## (e) Evaluator notes (documented, not exploited)

* τ_pass bookkeeping: a polar mirror that traps a few markers removes them from the passing set and raises τ_pass without any
  change in retention (cone −0.05: retained 413 = host 412 at 480 markers, score +3 %). Compare retained counts, not τ_pass.
* The wall/wire loss labels are not physics: a conductor placed in an escape path converts "wall" into "wire" one-for-one.
* `hopfmirror` clearance is set by the outer strand's crowding in the hole (0.016 for η 1.0 δη 0.06 at any N 6–24; 0.040 only
  for N 4 δη 0.2), so the global `wire_radius` couples an auxiliary torus to the main coil's buildability.
* `agent/composition.propose()`: 'add'/'replace' jitter ±0.12 sends donors at ρ ≥ 0.9 outside the ball (6/40 per round);
  'duplicate' flips the current sign at random and only translates along one axis.
* iota aliasing (brief §6) confirmed on the TF designs: true ι −4.6…−219 where `--surfaces` would read ≈ 0.
* The two 480-marker confirmations used the default (matched) seeds, not the held-out 777/4242/999 set — chosen deliberately so
  the cone-vs-host comparison is paired.

## (f) What I would try next

* Nothing more on plugs of the magnetosphere's null; the win is in the separatrix shape (fable-rings) and in reducing the
  near-coil seed losses, which is a seed-shell/geometry property (a smaller-radius stack keeps the same wire fraction).
* A symmetric-pair operator for the structural search (add/duplicate/mutate mirrored pairs and coaxial stacks) and a
  proxy that includes the retained-count anatomy; run it on the belt leader with sigma ~0.1 on the pair parameters only.
* If a toroidal field is really wanted for plasma reasons, the only geometry that does not cut the closed lines is an
  axial conductor with returns outside the wall — not available as a closed family (the `solenoid` leads break axisymmetry).
