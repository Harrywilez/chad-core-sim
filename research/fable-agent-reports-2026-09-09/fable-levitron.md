# fable-levitron — phase-2 report (2026-09-09)

Family: **Levitron** = nested `hopfmirror` toroidal-field pair (η ≈ 0.85–0.90, N 24, δη 0.06, wire at the 0.45 × clearance
cap) + one internal axis `circle` ring at ρ ≈ 0.40 carrying 2–3 × the strand current. A tokamak-like field whose "plasma
current" is the levitated ring: real rotational transform, buildability 1.0. Phase-2 work: 48 screening evaluations
(b2–b6), 2 confirmations at 480 markers / held-out seeds, 1 `--surfaces` run, 4 unaliased-ι runs (`tools/iota.py`), plus
two local diagnostics (marker loss locations, poloidal-contour maps) that explain what limits the family.
Notes with every hypothesis → result: `notes/fable-levitron.md`. Generator: `designs/fable-levitron/gen.py`.

## (a) Top 5 (fixed 120 T; 160 markers matched seeds; "@480" = 480 markers, held-out seeds 777/4242/999)

| # | design (file in `designs/fable-levitron/`) | score | score_plasma | τ_all / τ_pass | build | closed | S120 | unaliased ι (ring → nest edge) |
|---|---|---|---|---|---|---|---|---|
| 1 | `L_e0.90_N24_r0.400_I+2.00.json` (phase-1 leader) | 89.9 (**94.6 @480**, ref ver1) | 15.6 (`--surfaces`, ref surf1 — artifact, see (c)) | 85 / 95 (93 / 97 @480) | 1.00 | 0.44 | 0.31 | −1.6 (r 0.07) → −0.56 (0.12) → −0.21 (0.17); open beyond 0.21 |
| 2 | `L_e0.90_N24_r0.400_I+2.00_vf0.10at0.80z0.40.json` (+ weak opposing VF pair) | 91.3 (**92.9 @480**) | not run | 87 / 96 (94 / 92 @480) | 1.00 | 0.36 | 0.32 | −1.33 (0.07) → −0.39 (0.12); open beyond 0.16 |
| 3 | `L_e0.85_N24_r0.400_I+3.00.json` (η 0.85, ring 3 ×) | 94.7 (**90.3 @480**; 91.3 in the surf2 run) | 15.9 (`--surfaces` surf2 — same artifact) | 93 / 97 (88 / 92 @480) | 1.00 | 0.48 | 0.36 | −3.3 (0.07) → −1.0 (0.12) → −0.50 (0.17) → −0.26 (0.21); open beyond 0.26 |
| 4 | `L_e0.85_N24_r0.400_I+2.50.json` | 95.1 | not run | 94 / 96 | 1.00 | 0.46 | 0.36 | (same family as #3; not traced) |
| 5 | `L_e0.85_N24_r0.400_I+2.00.json` | 89.3 | not run | 90 / 89 | 1.00 | 0.40 | 0.33 | (not traced) |

Runners-up at the noise level: η 0.80 I 3.0 87.8; η 0.85 I 3.5 86.4 (τ_c 96.8, build 0.89); R 0.44 I 2.0 85.8; R 0.38 I 2.0 85.6;
R 0.40 I 1.5 85.1. **All of these are one plateau**: score 85–95, τ_c 85–97, buildability 1.00, S120 0.29–0.36. At 480 markers
the three confirmed designs are 90–95; nothing in this family is separable from the phase-1 leader by more than the ±5 %
held-out noise. The honest family number is **score ≈ 92 ± 4 (τ_c ≈ 92, build 1.0)**.

Loss channels (leader, 160 markers, evaluator counts: retained 50, wall 35, wire 75; local 40-T diagnostic with the same
seeds and field, final positions): the "wire" losses are **the inner strand shell (54 of 67), not the ring (6)** — hits
cluster at the top/bottom of the inner torus (poloidal angle 45–135° about the tube centre, |z| ≈ 0.31). Wall losses
(33, |z| median 0.71, ρ ≈ 0.4, lost in < 1.6 T) are markers seeded in the tube's hole and above/below the tube — on the
ring's vertical lines. By start distance from the ring: d < 0.1: 11/17 retained; 0.1–0.2: 39/56; 0.2–0.3: 10/35;
> 0.3: 0/52. I.e. the confined plasma is the nest of radius ≈ 0.19–0.22 around the ring (≈ 45 % of the seed shell);
everything else is lost in ≤ 12 T. Same picture for #3 (retained 0.36).

## (b) What I learned (confirmed / falsified)

1. **Unaliased transform (item 1)**: the nest has |ι| from 1.3–3.3 next to the ring to 0.2–0.3 at its edge — far above
   the 0.02 shorting threshold. The hybrid's ring-dominated twin (`a4_dip40_pair80_040`, ratio 2.5) has ι −4.0 → −0.31,
   same nest size 0.21: more ring current buys transform, not nest size. The evaluator's `--surfaces` ι_edge (0.0007 for
   the leader) is wrong for two stacked reasons: (i) modulo-1 aliasing of the inner values (−1.6 reads −0.09, −0.56 reads
   0.007), (ii) its "outermost surface" (r 0.37 from the axis) is a near-wall toroidal-solenoid line with ι ≈ 0 that is not
   part of the ring nest (iota.py: everything beyond 0.21 outboard is wall-bound). → score_plasma 15.6 is an artifact.
2. **(R, I) plateau (b2, 14 evals)**: R 0.36–0.44 × I 1.5–2.5 at η 0.9 gives 69–86, no ridge; I 3.0 at η 0.9 costs
   buildability (0.91). The predecessor's "ring ≈ 2 × strand" is confirmed as the buildability-limited optimum at η 0.9.
3. **η 0.85 (b5/b6)**: clearance 0.0235 (a 0.0105) lets the ring run at 2.5–3 × with build 1.0 → τ_c 93–97, the top of the
   plateau (95.1 / 94.7 screening, 90.3 @480). η 0.80 (a 0.0117, ratio 3–4) 84–88: no further gain. Weak (+5 %) and at the
   noise level, but consistent across 5 designs.
4. **Nest size is set by the pair's residual vertical field, not by the wall**: the bare hopfmirror pair carries a coherent
   B_z ≈ −0.045 × B_tor across the whole tube (cached grid of `agent2_r90N24`). Poloidal-contour maps: the ring's
   contours are circles (z_max = d, ρ_min = R − d) out to d ≈ 0.18; beyond that they open through the hole. Nest radius
   0.19 (ratio 2) → 0.22 (ratio 3). Closed fraction grows with ring current (0.28 → 0.53 from I 1.5 to 3.0) but τ does not.
5. **Scaled tube = FALSIFIED hard (b3, 11 evals, 16–43)**: shrinking the pair to sit inside the wall (scale 0.82/0.90) with
   the ring centred in the tube keeps the closed fraction (0.43–0.52) but retention collapses (S120 0.03–0.11). Reason:
   the "nested pair" is two different tori — inner strand ρ 0.12–0.82, z ±0.35; outer ρ 0.10–1.00, z ±0.45 — and the
   ∇B-drift excursion of nest-edge markers (Δ ≈ πρ_L/(2|ι|ξ) ≈ 0.03/|ι| ≈ 0.1–0.15) needs ≳ 0.15 of vertical clearance
   between the nest top and the inner-torus top: 0.16 in the leader, 0.10 when scaled. **Closed lines ≠ confined markers.**
6. **Vertical field (item 4, b4 + b6, 8 evals) = no gain either way**, exactly as the flux-function argument predicts for a
   ring inside a strand shell: an opposing VF (magnetosphere sign) closes the outboard contours *through the hole*, i.e.
   through the inner strand shell (closed fraction falls 0.44 → 0.36/0.28/0.19/0.18 for vf 0.1/0.2/0.3/0.4; score 91 → 76);
   a same-sense VF enlarges the nest (closed 0.49/0.51) but the added surfaces are drift-lost and the core with them
   (61.5 / 43.2). The magnetosphere trick needs a conductor-free hole; the Levitron does not have one.
7. Falsified small moves: δη 0.04/0.05 (raise the inner-torus top: 61–80, clearance loss wins); 2-turn vertical ring
   pair z ±0.05/0.08 (66 / 84); ring offset z ±0.06 (49 / 57 — symmetric penalty, the drift is two-sided); 2-turn ring
   stack (only matters for J, and the ring intercepts just ~4 % of markers — item 3 is moot).
8. Ceiling: the inner torus (hole 0.12, top 0.35) contains ≈ 75 % of the seed shell, the nest ≈ 45 %; even perfect nest
   retention gives τ ≈ 240 × 0.45 ≈ 110. The Levitron cannot reach the magnetosphere's 178 without a conductor-free
   hole, and a conductor-free hole means no toroidal field → no transform. The two leading families are complementary
   ends of the same trade-off.

## (c) Evaluator notes (documented, not exploited)

* `--surfaces` on ring-dominated fields: ι aliased mod 1 **and** the "edge" picked from near-wall ι ≈ 0 solenoid lines
  outside the nest → `shorted = false`, factor 0.17, score_plasma 15.6 for a design whose confined nest has |ι| 0.2–3.
  Quote both: evaluator ι_edge 0.0007 (leader) vs iota.py −0.21 at the last closed line. A fix would be to take the
  edge as the outermost surface *that encloses the axis found by iota.py's continuous θ*, and unwrap with the true ι.
  Second data point (surf2, `L_e0.85_N24_r0.400_I+3.00`): frac_surface 0.50, evaluator profile (r, ι) = (0.044, −0.34),
  (0.074, +0.09), (0.104, +0.05), (0.134, +0.35), (0.167, −0.33), (0.201, −0.20), then (0.397, −0.004), (0.403, −0.0008)
  → ι_edge −0.0008, factor 0.174, score_plasma 15.9. The first six values are exactly the iota.py profile
  (−3.3, ≈−2, −1.0, −0.65, −0.5, −0.26) taken mod 1; the last two are the near-wall solenoid lines at 0.40 from the axis
  (iota.py: open beyond 0.26). With ≥ 50 % surfaces and true ι_edge ≈ −0.2–0.26 this design would be exempt from the cap
  (score_plasma = score ≈ 90) if the edge were measured on the nest.
* Wire cap a = 0.45 × clearance makes the strand shells 90 % opaque (gap between wires 0.0023): the inner torus is an
  effective vessel wall. That is physics of the declared conductor, not a bug, but it means "wire" counts in
  `particles.counts` are strand hits, not ring hits, for this family.
* No ROI leak in any design here: roi_frac 0.57–0.59 and B_median 1.82–1.84 T throughout (the scaled tubes had
  roi_frac 0.38–0.54 with *lower* buildability — the opposite of a leak). Currents 104–114 kA per strand.
* Cross-machine / seed noise as advertised: 89.9 → 94.6, 91.3 → 92.9, 94.7 → 90.3 at 480 held-out markers.

## (d) What I would try next

* A **thinner inner shell**: the loss is the inner-torus top. A hopfmirror pair whose inner strand is *fatter* (top at
  z 0.45 like the outer one) is not available in this family — a `torushelix`-style toroidal-field winding on ONE torus
  (R0 0.56, r 0.44) with the ring inside would remove the inner shell entirely (tested family-wise by nobody yet).
* Ring ratio 2.5–3 at η 0.85 with a 2-turn ring stack (pitch 0.04) to hold build 1.0 at ratio 3.5–4 — τ_c 97 at I 3.5
  suggests another +3–5 %, i.e. still on the plateau.
* Confirm the plateau top with two more 480-marker runs (η 0.85 I 2.5, η 0.80 I 3.0) if the phase-3 question is "which
  file", not "which family".
