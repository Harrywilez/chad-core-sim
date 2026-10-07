# fable-weave — phase 2 report (2026-09-09)

Direction: the woven chiral meshed torus (`hopfmesh`) — stellarator surfaces, no conductor inside the plasma. Question: which
(N, δη, η, ε, periods) keeps ≥ 20 % nested surfaces AND buildability ≥ 0.9, and what is its best honest score_plasma?

Answer: **the same-strand-limited high-N family, N 40–48 at δη 0.016–0.018 (η 0.60, ε 0.70, n 4)** — buildability 0.93–1.0,
5/16 surfaces (0.31), ι −0.13 → −0.18, **confirmed score = score_plasma 41.4 at 480 held-out markers** (`c2_N48_de016`) — twice
the reference woven torus (20.0 / 14.6 held-out) and above the Levitron target (30). Every leader below is exempt from the plasma
cap (surfaces ≥ 0.25 and |ι_edge| ≥ 0.13 ≫ 0.02), so score_plasma = score.

## (a) Top 5 (fixed 120 T; *conf* = 480 markers, held-out seeds 777/4242/999, cloud, with `--surfaces`; screen = 160 m matched, vm)

| # | design file | score (conf) | score_plasma (conf) | screen 160 m | τ_all / τ_pass (conf) | build | closed | frac_surface | r_out | ι axis → edge |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `designs/fable-weave/c2_N48_de016.json` (N 48, δη 0.016, η 0.60, ε 0.70, n 4) | **41.4** | **41.4** | 45.9 | 38.0 / 51.6 | 0.93 | 0.35 | 0.3125 (evaluator) | 0.119 | −0.127 → −0.182 |
| 2 | `designs/fable-weave/c3_N24_de025_e65.json` (N 24, δη 0.025, η 0.65) | 36.4 | 36.4 | 51.5 | 32.2 / 41.2 | 1.00 | 0.19 | 0.25 (evaluator; one seed from the cap) | 0.110 | −0.100 → −0.139 |
| 3 | `designs/fable-weave/c2_N40_de018.json` (N 40, δη 0.018, η 0.60) | 31.1 | 31.1 | 44.6 | 29.4 / 32.8 | 1.00 | 0.49 | 0.3125 (evaluator) | 0.120 | −0.134 → −0.200 |
| 4 | `designs/fable-weave/c5_N24_de025_e65_vfp30.json` (N 24, δη 0.025, η 0.65 + VF ring pair ρ 0.55 z ±0.50, I +0.30) | **59.0 (screen only, unconfirmed)** | 59.0 | 59.0 | 51.1 / 68.1 | 1.00 | 0.31 | 0.3125 (local) | 0.128 | −0.075 → −0.129 |
| 5 | `designs/fable-weave/c1_N24_de020.json` (N 24, δη 0.020, η 0.60; the "big state") | 41.3 (160 m cloud, official `--surfaces`) | 41.3 | 38.9 (vm) | 47.4 / 64.4 | 0.75 | 0.53 | **0.4375 (evaluator)** | 0.149 | −0.100 → −0.218 |

Runners-up (screen only): `c4_N40_de020_e65.json` 48.5 (build 1.0, 4/16 local), `c5_N24_de025_e65_vfp60.json` 49.9 (6/16 local,
r_out 0.149, ι −0.06 → −0.13), `c1_N22_de025.json` 41.3 (4/16 local).

Also confirmed at 480 held-out by the ref agent: `b1_N20_de025` 27.3 (surfaces 6/16 local, big state). Expect any screen number of
this family to lose 10–30 % at 480 held-out markers (few markers land in the ρ 0.33–0.6 confined region).

## (b) What was learned

1. **Clearance = 0.66 δη for every N** (the inboard over/under separation), and the same-strand spacing ceiling is ≈ 0.50/N
   (0.0124 @N40, 0.0104 @N48, 0.0089 @N56). Current per turn ∝ 1/N, so buildability ∝ N δη² until the ceiling. The optimum is
   where 0.66 δη_min(N) meets 0.50/N: **N 40–48, δη 0.016–0.018** (build 0.93–1.0). τ_c does not fall with N (42–49 at N 36–48,
   the reference's 47 at full buildability); the ripple penalty of the predecessor's N ≤ 14 designs saturates by N ≈ 18.
2. **Surfaces vs N (the requested map).** The hypothesis "surface volume falls with ripple ∝ 1/N, crossover at N 14–16" is
   falsified: r_out is set by the weave state, not by N. Two states exist (Poincaré fans, `surf_local.py fan`):
   *big* — axis at R 0.485, nested surfaces spanning ρ 0.33–0.66 (6–7/16, r_out 0.125–0.149): N24 δη0.02, N20 δη0.025, N18
   δη0.025, N16 δη0.03 (5/16); *small* — axis at R 0.39–0.42, surfaces ρ 0.33–0.53 (4–5/16, r_out 0.10–0.12): N22/N24/N28 δη0.025,
   N28/N32 δη0.02, all N ≥ 36. The zone fraction of the zone weave predicts the state (big ⇔ zones ≤ 0.93, small ⇔ ≥ 0.97;
   zones ≈ g(N δη), threshold N δη ≈ 0.5). In the small state r_out grows as the shell separation shrinks (sep 0.0165 → 4/16,
   0.011 → 5/16). The 20 % cap fails only at N10 δη0.04 (3/16), η ≥ 0.65 at N ≤ 22 δη 0.025, n5/n6, ε 0.55 (2–3/16).
3. **Fatness η 0.65–0.70** raises the raw score 30–40 % (τ_pass 52–65; more of the seed shell inside the tube) but pulls the axis
   inboard (0.41 → 0.36 → 0.31, the inboard conductor is at 0.28/0.25/0.22) and shrinks the surfaces (r_out 0.11 → 0.06):
   η 0.65 passes only at 4/16 (N24 δη0.025, N40 δη0.02, N44 δη0.019); η 0.70 fails (2/16). The raw gain is tube confinement,
   not surfaces. Hypothesis 3 (η 0.65–0.75 at the best point) is falsified as a surfaces lever; η 0.65 is a legitimate but
   marginal raw-score lever.
4. **Modulation** (hypothesis 2) falsified for the nested state: at N24 δη0.025, n3 4/16 (ι −0.25 → −0.30), n5 3/16, n6 2/16,
   ε 0.85 4/16 (ι doubles to −0.21 → −0.28, score 33), ε 0.55 3/16; scores 33–43 vs 36.7 at n4 ε 0.70. RUNGS §0's "n5 = 50 %"
   holds only for the woven N24 torus. `helical_mode` poloidal was already falsified in phase 1.
5. **A vertical field is the axis lever — and it raises τ too** (new): a ring pair at ρ 0.55, z ±0.50 (r 0.74, outside the ROI
   shell; roi_fraction and the strand current are unchanged, so no leak) with +0.3 × strand current moves the axis of the η 0.65
   leader from 0.359 to 0.387, the surfaces from 4/16 to 5/16 (r_out 0.110 → 0.128) and the screen score from 51.5 to **59.0**
   (τ 51/68, build 1.0); +0.6 gives 6/16 (axis 0.430, r_out 0.149, ι −0.06 → −0.13) and 49.9; −0.3 gives 2/16 and 40.8.
   Confirmation at 480 held-out markers is the first thing to do next (expect ≈ 0.7–0.9 × 59).
6. Held-out 480-marker scores are 10–30 % below the matched 160-marker screen for the whole family (ref 20.0 → 14.6; b1_N20_de025
   36.3 → 27.3; c3_N24_de025_e65 51.5 → 36.4; c2_N40_de018 44.6 → 31.1; c2_N48_de016 45.9 → 41.4). Rank by confirmations only.

## (c) Evaluator notes / artifacts (documented, not exploited)

* `--surfaces` uses 16 seeds on the outboard ray from the axis; frac_surface is quantised in 1/16 and the 20 % threshold is
  effectively 4/16 = 0.25. Several good designs sit exactly at 4/16 (one seed from τ_pol ≈ 4–5 transits, i.e. a ×0.1 cliff).
  A finer seed count or an area-weighted measure would make the cap far less noisy.
* The axis finder (`find_axis`) starts its fan at R0 − 0.55 r_tube = 0.44 and walks inboard by ≤ 0.03 per refinement; for the
  fat η 0.70 tube the true axis (0.31) is at the edge of what it can reach. It found the right axis in every case I checked with
  a full midplane fan, but the "reach" (axis → outboard strands) then spans 0.55–0.66, so the 16 seeds sample a region that is
  75–85 % open by construction.
* The surfaces trace uses wall = 1.05 and wire_clear 0.02, the orbits the vessel r = 0.82: a surface can be counted although it
  crosses the wall (the big-state surfaces reach ρ 0.66 < 0.82, so it did not matter here).
* Local reproduction: `designs/fable-weave/surf_local.py surf` runs the evaluator's own routine and matched the queue's
  frac_surface / r_out / ι exactly on the four designs with official runs. I used it (one process at a time, niced) as a
  pre-screen; ~25 local traces, 15–200 s each. No ROI leak in any design (roi 0.65–0.67, B_median 2.2 T throughout).

## (d) What to try next

1. Confirm `c5_N24_de025_e65_vfp30` (59.0 screen) at 480 held-out with `--surfaces`; then VF +0.3…+0.6 on `c2_N48_de016` and
   `c4_N40_de020_e65` (build 1.0, raw 48.5), and a finer VF scan (0.2–0.5, ring radius 0.45–0.65, z 0.45–0.6) — if the surfaces go
   to 6–7/16 with τ intact, 45–50 confirmed is plausible.
2. A big-state design at full buildability: the zone threshold N δη ≈ 0.5 gives N16 δη0.03 (build 1.0, 33.3), N18 δη0.0275,
   N20 δη0.025 (0.95); add VF to push the axis further out; try η 0.62–0.63 (compromise between tube volume and axis position).
3. Understand the state switch: the free (ramp) stretches are outboard; a deliberate "woven outboard, nested inboard" pattern
   (e.g. δη varying with the poloidal angle, if the generator allowed it) could give the big state at any N.
4. 480-marker confirmation of `c1_N24_de020` (τ_c 52, the largest surfaces) and of `c4_N40_de020_e65`.

Budget used: 41 screening evaluations (c1 13, c2 8, c3 9, c4 8, c5 3) + 4 `--surfaces` (3 combined with the 480-marker
confirmations, 1 on c1_N24_de020) + 3 confirmations; ~70 minutes wall-clock. Notes: `notes/fable-weave.md` (per batch), `notes/fable-weave_surf.jsonl` (local
surfaces), `notes/fable-weave_quick.jsonl`. Generators: `designs/fable-weave/gen_p2.py`, `mk.py`; diagnostics `surf_local.py`.
