# fable-cnt report (phase 2, 2026-09-09, 04:55–06:20 UTC) — stellarator surfaces from a few planar circular coils

Direction: CNT-type interlocked coils + PF pair, and heliac-type sets (central ring + planar TF coils on a helical axis).
Budget used: 45 screening evaluations (160 markers, 120 T), 2 confirmations at 480 markers, 2 `--surfaces` jobs
(both finished — see §3), ≈ 40 quick.py checks and one clearance-only scan. Designs in
`designs/fable-cnt/` (generator `gen.py`, scan driver `scan.py`, closed-line map `diag.py`); running notes in
`notes/fable-cnt.md` (a clock-correction note explains the fast section timestamps); quick-check log
`notes/fable-cnt_quick.jsonl`. Result JSONs: `queue/done/fable-cnt__*`.

## 1. Top 5 (fixed 120 T, matched seeds; build = 1.00 for all; `score_plasma` where a surfaces job ran)

| # | design (`designs/fable-cnt/…`) | score | score_plasma | τ_all / τ_pass | build | closed | S120 | notes |
|---|---|---|---|---|---|---|---|---|
| 1 | `b6_hax_Rc54_N30_rk24_M3_d12.json` | **54.9** (54.1 with `--surfaces`, cloud) | **54.1** (factor 1.0) | 40.6 / 74.2 | 1.00 | 0.33 | 0.14 | 30 planar circles r 0.24 on the helix ρ = 0.54 + 0.12 cos 3φ, z = 0.12 sin 3φ; clearance 0.0251 = chord cap; 334 kA |
| 2 | `b5_hax_Rc52_N30_rk24_M3_d12.json` | 48.6 (**41.2 @ 480 markers**) | — | 38.1 / 62.0 (33.9 / 49.9 @480) | 1.00 | 0.41 | 0.13 | same with Rc 0.52 |
| 3 | `b5_hax_Rc54_N32_rk24_M3_d12.json` | 44.8 | — | 33.9 / 59.2 | 1.00 | 0.41 | 0.12 | N 32 |
| 4 | `b6_hax_Rc55_N30_rk26_M3_d12.json` | 42.9 | — | 34.1 / 54.0 | 1.00 | 0.32 | 0.12 | fatter coils r 0.26 (clr 0.0247) |
| 5 | `b5_hax_Rc52_N32_rk24_M3_d12.json` | 42.8 | — | 35.2 / 52.1 | 1.00 | 0.36 | 0.12 | |
| ref | `b3_hax_Rc50_N32_rk22_M3_d12.json` (baseline) | 37.2 / 38.0 (34.2 @480) | **38.0** (factor 1.0) | 31 / 44 | 1.00 | 0.45 | 0.10 | 100 % nested surfaces, ι 0.36 → 0.42, r_out 0.19 |

Best CNT-type design (for the record): `b2_cnt_a50_d10_t85_pf95z30i-70.json` — 23.9 (τ 34.8 / 16.5, build 1.00, closed
0.12, S120 0.07): IL coils a 0.50 at (∓0.10, 0, 0) tilted ±42.5° about x, PF pair R 0.95 z ±0.30 at −0.70 × I_IL.

The winner is a **helical-axis stellarator made of 30 identical planar circular coils** (a heliac without its central
conductor): no conductor inside the plasma, no weave, buildability 1.00 at the declared wire (0.45 × the lone-circle
chord cap; the true coil–coil clearance equals the cap), and it is plasma-legit: the baseline's `--surfaces` run gives
frac_surface 1.00, ι_axis 0.355 → ι_edge 0.422 (no aliasing: |ι| < 0.5), shorted, `surfaces_ok`, score_plasma = score.
With 160-marker noise (c1/c2 show the family reads ≈ 10–15 % high at 160), the honest number for #1 is ≈ 47–55; it
beats the woven chiral torus (36.3 → 27.3 @480; reference 20.0) and the Levitron reference with surfaces (29.6) as the
best plasma-legit stellarator seen in this search, but it is far below the poloidal-field families (magnetosphere 178,
Levitron 90, whose score_plasma is capped or aliased).

## 2. What was learned (confirmed / falsified)

1. **CNT-type (2 interlocked circles + PF pair)**: closed lines exist only with a *negative* PF current (opposing the IL
   pair's net vertical dipole 2Iπa² cos θ/2), −0.5 to −0.8 × I_IL, PF coils large and low (R 0.90–0.95, z ±0.30–0.40).
   Halving the IL offset d 0.20 → 0.10 tripled the score (7.7 → 23.9); tilt 70–90° is flat within noise; a 0.50 gives
   build 1.00. **Negative**: it cannot become competitive — the CNT plasma is a thin twisted torus (≈ 10 % of the seed
   shell) with low ι; S120 ≤ 0.07 for every variant.
2. The predecessor's "true CNT frame" (`cnt2`, coils' common axis along z) with z-coaxial PF coils gives 0.00 closed
   lines for any PF current: the plasma torus circles the *bisector of the two coil normals* (the y axis there), so those
   PF coils were perpendicular to it. The `cnt` frame (common axis x, bisector z, PF along z) is the physical CNT.
3. **Heliac with a central ring: negative.** With 8–16 TF coils (≤ 5 per period, tilts up to 42°) no line closes; with
   24–32 coils and a ring at ±3 × I_TF the closed fraction is 0.00–0.04 — the ring destroys the surfaces the TF coils
   make by themselves. Removing the ring gives closed 0.40 → the helical-axis stellarator.
4. **Helical-axis optimum**: M 3, δ 0.12 (δ 0.08/0.10/0.16 all worse; M 2 worse; M 4 collides), N 30–32 (24 → more
   current per coil, build 0.78; 36–40 no better), coil radius as large as the clearance allows: rk 0.22 → 0.24 raised the
   score 37 → 43–55 once Rc was moved out to 0.52–0.54 so that neighbouring tilted coils clear each other (clearance =
   chord cap). Rc 0.56 collapses (30.5): the wall r 0.82 cuts the outboard tube (ρ_out 0.92). Rc 0.45 (tube fills more of
   the seed shell, closed 0.56) is worse (24): weaker inboard field (B_med 0.39 T) and higher current.
   **Coil–coil clearance, not the field, limits this family** (rk 0.26 at Rc 0.50–0.52: τ_pass 50–60 but build 0.2–0.8).
5. A ±0.10 vertical field (PF pair R 0.90, z ±0.40) changes nothing (37.6 / 36.2 vs 37.2).
6. τ_pass > τ_all in every helical-axis design (surfaces hold the passing particles; ripple-trapped ones are lost,
   R_mirror 2–3 with 30 coils) — the opposite signature to the mirror families.

## 3. Does the `--surfaces` analysis work for this family? — Yes

Baseline job `fable-cnt__b3s__b3_hax_Rc50_N32_rk22_M3_d12__199eaf95`: axis found at R 0.595 (the helix at the section's
φ), frac_surface 1.00, frac_island_chaotic 0, outermost surface r 0.189, ι_axis 0.355 → ι_edge 0.422 (|ι| < 0.5, so no
rotation-number aliasing; consistent with the M 3 helical excursion δ/Rc ≈ 0.24), shorted = True, surfaces_ok = True,
plasma factor 1.0, **score_plasma 38.0 = score**. The #1 design (`fable-cnt__b6s__b6_hax_Rc54_N30_rk24_M3_d12__01f895c3`) confirms it: axis R 0.633,
frac_surface 1.00, islands/chaotic 0, outermost surface r 0.209, ι_axis 0.305 → ι_edge 0.369, shorted, surfaces_ok,
**score_plasma 54.1 = score** (the two surfaces jobs took 8 and 30 min on the cloud; 22 000–23 000 conductor points). `tools/iota.py`
on the baseline exceeded 300 s on the loaded machine and was abandoned; the evaluator's own ι is unaliased here anyway.

## 4. Evaluator notes (documented, not exploited)

* The 25³ field grid (dx 0.0725) cannot resolve a ring's 1/d field within ~0.1 of it: at 0.03 above the heliac ring the
  interpolated field was 97 % toroidal where the ring alone should give 1.6 × the TF field. Every "internal ring" family
  (Levitron, heliac) is evaluated in a smoothed version of its own field.
* Lone-circle chord cap: the helical-axis clearance 0.0230–0.0251 (rk 0.22–0.24) and the CNT IL clearance 0.0523 (a 0.50)
  are the 0.105 × radius artifact; for the top designs the true coil–coil distance equals it (checked with `validate` at
  rk 0.26, where it drops below the cap). wire_radius was declared at 0.45 × the cap, as instructed.
* 160-marker noise is large for this family: identical-physics designs (baseline ± PF 0.1) read closed 0.25 / 0.45 / 0.59
  and the two 480-marker confirmations came in 8 % and 15 % below the 160-marker screens.
* `roi_fraction` 0.61–0.68 for 30–36 coils comes from the 0.08 conductor exclusion around many coils, not from a parked
  auxiliary; B_median_shell 0.5–1.1 T reflects that the inner third of the seed shell is outside the tube.
* `--surfaces` runtime 8–30 min for 30-coil designs makes it unusable for screening this family (≤ 1 per session).

## 5. What to try next

* The clearance frontier: rk 0.26–0.28 with Rc 0.54–0.55 and δ 0.10–0.12, N 28–30, checked first with `validate`
  (clearance ≥ chord cap keeps build 1.00); Rc 0.55 N30 rk 0.26 δ 0.10 (clr 0.0257) was written but not evaluated.
* Elliptical or D-shaped coils are outside the "planar circles" brief, but a second ring of 30 smaller circles (r 0.12)
  nested inside the tube at 2 × current could raise the inboard field without collisions — not planar-only-outside-plasma.
* M 3 with a small helical *pitch modulation* (coils bunched near the inboard crossings) to lower the ripple that loses
  the trapped particles (R_mirror 2.7 → τ_all ≈ 0.55 τ_pass).
* Confirm #1 at 480 markers on held-out seeds before quoting 55 (its surfaces result is in: score_plasma 54.1).
