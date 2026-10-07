# fable-rings — phase-2 report (2026-09-09): the magnetosphere, pushed to its floor

Direction: a multi-turn levitated ring inside an opposing external field (dipole + "anti-Helmholtz" = spherical
separatrix, every field line closed). Phase 1 ended at 178.5 (addendum) — its last batch actually finished after the
cut-off at **193–197** (R 0.45 n 2 / R 0.50 n 2 + Helmholtz pair, build 1.00). Phase 2: 41 screening evaluations,
3 confirmations (480 markers, held-out seeds 777/4242/999, grid 41), 2 `--surfaces`. Notes: `notes/fable-rings.md`;
generator: `designs/fable-rings/_gen2.py` (batches c1–c4).

## (a) Top 5 (fixed 120 transits, 160 markers, matched seeds; confirmation column = 480 markers, held-out seeds, grid 41)

| # | design (file in `designs/fable-rings/`) | score | τ_all / τ_pass | build | closed | S120 | wall / wire (of 160) | confirmed @480 m, grid 41 |
|---|---|---|---|---|---|---|---|---|
| 1 | **`c3_R40tri_belt_i92.json`** — R 0.40 3-turn triangle bundle (pitch 0.041, a 0.0184) + belt pair (ρ 0.90, z ±0.25) I −0.92/turn | **216.7** | 216.7 / 216.7 | 1.00 (J/J_REBCO 0.79) | 1.00 | 0.90 | 7 / 9 | **222.4** (S120 0.92; 441 held / 21 wall / 18 wire of 480) |
| 2 | `c2_R40tri_belt_i80.json` — same bundle, belt I −0.80 | 214.2 (216.8 at grid 41) | 214.2 / 214.2 | 1.00 | 1.00 | 0.89 | 9 / 9 | 220.1 (S120 0.91; 437 / 24 / 19) |
| 3 | `c3_R40tri_belt_i86.json` — same bundle, belt I −0.86 | 214.6 | 214.6 / 214.6 | 1.00 | 1.00 | 0.89 | 7 / 11 | — |
| 4 | `c2_R45n2_belt_i63.json` — R 0.45 2-turn column (pitch 0.046, a 0.0207) + belt I −0.63 | 210.8 | 210.8 / 210.8 | 1.00 (0.87) | 1.00 | 0.87 | 9 / 12 | — |
| 5 | `c1_belt_R45n2_p90z25_i60.json` — R 0.45 2-turn column + belt I −0.60 (the 4-circle "canonical" design) | 210.3 (cloud repeat 208.2; grid 41 209.5) | 210.3 / 210.3 | 1.00 | 1.00 | 0.86 | 8 / 14 | 216.3 (S120 0.89; 427 / 29 / 24) |

All conductors: `circle` family only; the belt coils sit at r = 0.934 (outside the wall 0.82, inside the ball); wire radius
= 0.45 × clearance (the bundle pitch, itself at the lone-circle chord cap 0.105 R); roi_fraction 0.92–0.93 and
B_median 1.5–1.6 T for every design (no ROI-mask leak: the belt is 0.2 outside the ROI). Also on the plateau (204–211):
R 0.35 2×2 bundle, R 0.38 2×2, R 0.40 3-turn column, R 0.45/0.38 triangles, 4-turn diamond, 3-coil belt.
score_plasma: `--surfaces` gave 208.2 (= score) for #5 and 14.1 (capped) for #1 with identical topology — an artifact,
see (d); only the raw score is claimed. Budget used: 41 screening + 3 confirmations (480 m) + 2 `--surfaces`, ~75 min.

## (b) The loss-channel finding — where the last 12–16 % go

Reconstructing the seed pool for each finished job (positions are deterministic) and computing ψ = ρA_φ at every seed:

1. **Wall losses = the polar column.** For the phase-1 leader (22 wall of 160) *every* wall-lost marker had
   ψ/ψ_max < 0.022 — seeds within ρ ≲ 0.35 of the axis at |z| 0.2–0.5, where |B| is 0.1–0.35 of the median — and only 3 of
   125 retained markers sat in that flux band. Pitch angle is irrelevant (|ξ₀| of lost = retained). Orbit-free tracing of the
   axial null and the equatorial ψ = 0 crossing showed why: **the real separatrix is oblate** (a ring's on-axis field is weaker
   than the point dipole's, its equatorial far field stronger). With the Helmholtz pair (0.80, ±0.40) the equator touched the
   wall (ρ_s 0.81) while the polar null sat only 0.05 above the seed shell (z 0.63); a weaker pair opened the equator first
   (closed 0.78), a stronger one pulled the null into the shell. The pair couples both.
2. **Fix = a belt pair** (ρ 0.90, z ±0.25): more opposing field at the outboard equator, less on the polar axis → a *spherical*
   separatrix at r ≈ 0.71 (z_null 0.71 = ρ_s 0.72), half-way between shell (0.58) and wall (0.82); weak-field overlap with
   the shell 4 % → 0. Result: wall losses 22 → 7–9, all nine belt designs beat all five Helmholtz ones (+10–15 points).
3. **What remains is two floors.** (i) Polar column, 5–6 %: the same ~9 seeds (ρ ≤ 0.15, |z| 0.37–0.52, ψ < 0.015 ψ_max) die in
   every design within 5 transits. Near the axis ψ = ρ²B_z/2 (≈ 0.0013 T m² at ρ 0.08) is *below* the canonical-momentum
   excursion m v ρ/q (≈ 0.002 T m²), so p_φ conservation cannot hold them and μ fails at the null; holding ρ = 0.1 would need
   B_z(axis) > 0.5 T in the hole, above the 0.4–0.6 T any ring gives at the 2 T-rms normalisation (smaller rings, bundles and
   tall columns changed Bax/median by ≤ 20 %: no measurable gain). (ii) Conductor, 6–8 %: seeds 0.04–0.09 from the coil.
   Grid 41 (same seeds) leaves the bundles' hits unchanged (9 → 9, 14 → 11) but halves a lone ring's (15 → 8, τ 199 → 214):
   for multi-turn coils the channel is real — lines encircling one turn thread the gap between turns (0.005, less than a wire
   diameter) — plus ~5 % of seeds effectively born at the conductor (the 0.03 seed filter uses the interpolated distance map,
   which over-estimates distance near a wire). A thinner declared wire (a × 0.93) changed nothing.
   Floors: 240 × (1 − 0.06 − 0.07) ≈ 209 at 160 markers; the confirmed 222 is the same physics with the held-out pool.

## (c) Lessons (confirmed / falsified)

1. CONFIRMED (H8): belt pair ≫ Helmholtz pair for a real ring — shape the separatrix, not just its size; the optimum is round.
2. CONFIRMED: the vertical-field optimum is a sharp threshold on equatorial closure (Helmholtz i51 156 → i57 181 → i63 193), then
   a plateau; for the belt the optimum is I_v ≈ 0.90–0.92 × I_turn for the R 0.40 triangle (i1.00 208, i1.08 204 fall off).
3. CONFIRMED: ring z-offset is strongly bad (+0.05 → 165.7, both channels doubled): z = 0 is not a free parameter.
4. CONFIRMED (predicted by the orbit-free proxy before spending evaluations): a polar pair (0.75, ±0.50) cannot close the
   equator (143, closed 0.58); the 1-s null/ψ-crossing diagnostic (`scratch proxy.py`, logic in the notes) ranked c1 correctly.
5. FALSIFIED (H9): conductor hits do not scale with the coil's near-zone volume (R 0.35 2×2 bundle 14–16 vs R 0.45 column 12).
6. FALSIFIED (H1 as posed): once the belt is there the coil shape (R 0.35–0.45, column/triangle/2×2/diamond, n 2–4) is a
   205–217 plateau; the compact triangle bundle is best by ~5 (fewer gap lines than a column: 9 vs 13–15 hits; a tall n 4 column
   is worst, 17–22 hits). Build 1.0 is free at R ≥ 0.40 with ≥ 3 turns (J/J_REBCO 0.6–0.8) — buildability was not the lever.
7. NULL (H10): thin wire (J → 1.0) does not reduce conductor hits; keep a = 0.45 × clearance.
8. Robustness: the leader confirms *upward* (222.4 at 480 held-out markers, grid 41; 220.1 and 216.3 for the runners-up);
   cross-machine 1 % (210.3 vm vs 208.2 cloud); grid 25 → 41 ±2 points.

## (d) Evaluator notes and artifacts (documented, not exploited)

1. **The plasma cap was *not* applied to this family, for the wrong reason.** `--surfaces` on the c1 winner returned
   score_plasma = score = 208.2 (factor 1.0): the Poincaré classifier reported frac_surface 0.56, outermost r 0.41, ι_edge
   −0.039 > 0.02 ("shorted"). `tools/iota.py --grid 41` on the same design: **"poloidal-only", Δφ/2π = 0.00** for every line
   (Δθ/2π −87…−193). B_φ of coaxial circles is identically zero; the ι is grid noise (a Cartesian grid breaks axisymmetry), and
   its sign flips along the profile (+0.009, −0.076, −0.129, +0.007 …). Addendum item 7 expected the cap; instead noise
   cleared the shorting threshold. The second run, on the confirmed leader `c3_R40tri_belt_i92` (identical topology), got
   ι_edge = +0.0034 → "not shorted" → factor 0.066 → **score_plasma 14.1** (τ_pol 14 transits) at score 215.2. Same physics,
   opposite verdicts, decided by the noise sign of a spurious ι. I quote the raw score and do **not** claim score_plasma;
   the physical argument is below.
2. Physics of item 7, for the record: in an axisymmetric poloidal field the ∇B and curvature drifts are toroidal — in the
   symmetry direction — so ions and electrons circulate on closed drift shells (third invariant) and no charge separation
   builds anywhere (∂/∂φ ≡ 0): there is no polarization field and no τ_pol loss. Planetary magnetospheres and LDX hold plasma for
   seconds. The cap is a model limitation for dipole-type fields; the transform argument does not apply to them. **Real
   caveats** the score does not contain: (i) the ring is a levitated superconducting coil *inside* the plasma — it needs a
   levitation coil, a cryostat in a 15 keV plasma and it is a particle sink (our 6–8 % conductor channel is exactly that);
   (ii) the interchange (flute) pressure-profile limit, p ∝ V^(−5/3): the profile may be no steeper — a transport, not a
   polarization, limit; (iii) the polar X-points (z ±0.68–0.71) scatter pitch angle and feed the polar column to the wall — a
   divertor region in a real device; (iv) the 2 T-rms vacuum normalisation ignores the high-β diamagnetism a dipole plasma
   reaches; (v) the belt coils (r 0.934) are outside the wall but inside the ball, as required.
3. Seed filter + coarse grid: seeds are kept if the *interpolated* distance map says > 0.03 from a conductor; near a wire
   the trilinear interpolation of a cone over-estimates distance, so 2–4 % of markers start inside/at the conductor (t < 0.05
   transits). Lone-ring conductor hits are halved by grid 41 (interpolation); bundle hits are not (real gap threading).
4. Grid 25 vs 41 moves the score by ±2 and nothing else; 160 markers give ±5 on this plateau (binomial on ~20 losses).
5. Stale leaderboard: the addendum's 178.5 was superseded by the predecessor's own last batch (193–197) that finished after
   the cut-off; the "b3 stragglers" were re-queued by the orchestrator and finished at the start of phase 2.

## (e) What next

* The ceiling of this family at this normalisation is ~225: two floors of ~6 % each. The polar column could only be reduced by
  a higher on-axis field in the hole, which the rms normalisation forbids — unless a legitimate design puts flux there
  (a same-sense inner coil on the axis is not axisymmetric-closed with `circle`; a small coaxial same-sense ring creates an
  X-line inside the shell, phase 1's double dipoles).
* The gap channel is physical for bare conductors but a real coil is a solid winding pack in a cryostat; a "porous" widely-spaced
  bundle would reduce it in the evaluator only — I did not pursue it (not honest engineering).
* If the grading moves to score_plasma with the cap fixed for poloidal fields, this family is the only plasma-legitimate
  one above 200; it needs the evaluator to recognise "poloidal-only" lines (as `tools/iota.py` does) instead of thresholding ι.
