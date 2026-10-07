# Composition screens on the corrected meshed torus — fable's runs (numba runtime)

Status: screening + partial validation. Evaluator versions are recorded per row in
`results/composite/ledger_fable_room_session.jsonl`; source hashes in `results/composite/source_hashes.json`.
All runs: reactor preset, B_rms(ROI) normalised to 2 T, 15 keV D⁺, seeds 12345/271828/31415
unless stated. My runtime is a separate machine with numba (`orbit_backend: numba`); it agrees
with the room's NumPy runtime to 0.3 % on the frozen cross-check (event 80 vs
`results/composite/crosscheck_index0_32m40T.json`: raw 7.7195 vs 7.7454, τ_c 18.25 vs 18.31,
closed 0.5833 both).

## 1. What was asked and what was run

Goal: let agents combine different components, then test many systems for confinement.
The combination mechanism exists in `ccsim.design` (component lists with scale / reflect /
rotate / translate / relative current) and mind-1's `agent/optimize.py` now proposes structural
moves. I ran physics-motivated compositions around the only design with verified flux surfaces,
the woven meshed torus II (`designs/composite/woven_m24_de015.json`: η 0.60, ε 0.70, n 4,
24 circuits, δη 0.015, conductor 0.0045 — score 17.2 at 480 markers, 44 % nested surfaces,
ι −0.10 → −0.20), plus a new `torushelix` family (classical l = 2 stellarator windings).

| composite (160 markers, repro1 unless noted) | τ_c | score | closed | note |
|---|---|---|---|---|
| woven torus, bare | 40.8 (480 m) | 17.2 | 0.58 | reference; 120 T |
| + axis ring R 0.64 in the tube, I = +0.1 / +0.3 / −0.3 / +1.0 | 23.9 / 14.3 / 18.8 / 7.5 | 10.0 / 6.0 / 7.9 / 3.2 | 0.54 / 0.35 / 0.50 / 0.06 | intercepts orbits, opens lines — monotonically bad |
| + vertical-field rings z ±0.46, R 0.62, I = 0.1 / 0.3 | 37.5 / 42.6 | 15.8 / 18.0 | 0.66 / 0.47 | = baseline within noise |
| + tornado cone plugs, same sense, I 0.3 / 1.0 (480 m) | 18.6 / 19.0 | 9.1 / 9.4 | 0.50 / 0.36 | halves confinement (AGENT.md §9) |
| manifest (mind-2, event 56): cones ±0.03 / ±0.10 | 40.7–44.2 | 19.8–21.5 | 0.48–0.64 | τ unchanged; score up only through the ROI leak (§3A) |
| manifest: cones +0.30 | 30.8 | 15.0 | 0.50 | negative control confirmed |
| manifest: rings R 0.3 at z ±0.5, +0.03 / −0.03 / +0.10 / −0.10 / +0.30 | 19.0 / 37.6 / … | 8.4 / 16.7 / … | 0.62 / 0.62 | +0.03 stopped at 40 T, −0.03 ran to 120 T: the ladder artifact (§3B), not physics |
| manifest: cones + rings +0.03 (three families) | see `manifest_results_repro1.json` index 11 | | | |
| ×0.85 woven torus alone | 19.3 (40 T) | 6.0 | 0.46 | scaled host for the helical windings |
| ×0.85 torus + l=2 helical windings n = ±4, I 0.3–1.0 | 43.9–51.0 (120 T) | 14.2–16.4 | 0.32–0.56 | see §2 — apparent gain is the ladder artifact |
| ×0.85 torus + windings n = 3 / n = 5 (2 conductors) | 16.6 / 37.3 | 5.3 / 12.0 | 0.44 / 0.52 | |

Full RESULT records for the 12 manifest cases: `results/composite/manifest_results_repro1.json`.

## 2. The helical-winding "doubling" is not confinement

At the default run ladder the ×0.85 torus stopped at 40 transits (τ_c 19.3) while every
n = ±4 composite ran to 120 transits (τ_c 44–51). Three checks say the windings do nothing
physically significant at these currents:

* Poincaré (`--surfaces`, v2.3): bare ×0.85 torus 44 % surfaces, r_out 0.124, ι −0.099 → −0.200;
  with windings (n 4, I 1.0) 44 %, r_out 0.134, ι −0.099 → −0.229. Same surfaces, same transform
  to within the seed resolution — an l = 2 winding carrying ~1 strand current against 46
  poloidal turns is a few-percent perturbation, and both handednesses gave the same τ.
* Fixed 40-transit horizon, same seeds, both designs: τ_all 16.0 vs 16.7, S40 0.13 vs 0.14,
  survivors 21/82 vs 23/85 — all of them markers seeded **inside the tube**.
* The 78 of 160 markers seeded in the hole, above/below the tube, or in the gap between tube
  and windings are lost within ~4 transits in both designs (S20 ≤ 0.06).

So the factor of two came from the 40/80/120 ladder: extension is keyed on S_T ≥ 0.10 and a
tail-share test, and designs sitting at S40 ≈ 0.10 land on either side of it. Confirmation runs
at 480 markers with a fixed 120-transit horizon and held-out seeds (777/4242/999) are in
`results/composite/comp_confirm*.log` as they complete.

## 3. Evaluator findings from these runs (beyond mind-2's radius fix and the fallback fix)

**A. ROI normalisation leak (open).** B_rms for the 2 T normalisation averages grid points with
r ∈ [0.20, 0.65] that are ≥ 0.08 from any conductor. A weak auxiliary parked in the weak-field
hole removes low-|B| points from the average; B_rms rises, the current for 2 T falls, J falls,
buildability rises 0.42 → 0.49 and the score 17.5 → 20.3 while τ_c is unchanged (manifest cones
±0.03). Any "weak auxiliary helps" result must be checked against this. v2.3 records
`roi_fraction`, `roi_points` and a mask-free `B_median_shell_T` so the leak is visible; a fix
(design-independent normalisation point set, or a normalisation over the confinement region
rather than the whole shell) needs a decision because it re-scales every historical score.

**B. Run-ladder threshold (open).** See §2. For any A-vs-B comparison run both at a fixed
horizon (`--transits 120`). A smooth alternative (extend on a fitted tail, or always run 120
and report S at 40/80/120) would remove the discontinuity at the cost of ~3× runtime for
short-lived designs.

**C. Current normalisation double-count (fixed in v2.3).** `ScaledGridField` normalises the grid
by the first winding's relative current, and I_max multiplied that current by rel_max again
(events 122/126/139). Now I_max = (current / |rel_first|) × rel_max; `MultiWinding.with_current`
preserves relative amplitudes (Joule power was silently equalised before); zero-current
components are disabled at build. `tests/test_evaluator_normalization.py`: uniform rescaling
of all relative currents and reordering of components leave I_max, buildability and the field
identical (one chaotic long-lived marker may differ at the 1e-3 level).

**D. Order sensitivity of the wire-distance grid (fixed in v2.3).** `FieldGrid` used every 2nd
conductor point for the distance grid; the ROI mask and the seed rejection then depended on the
concatenation order of the windings (2.5 % in I_max when two components were swapped). Every
point is used now (KD-tree); all cached grids must be rebuilt (`cache/designs/grid_*_n25.npz`).

**E. Half the markers start outside the plasma (open).** The seed shell r ∈ [0.18, 0.58] is
design-independent by intent (no design chooses its own test volume), but for a torus of tube
radius 0.31 centred at ρ 0.54 it puts ~50 % of the markers where there is no confining field.
Absolute τ values for toroidal designs are therefore dominated by this mismatch; ranking within
the toroidal family is still meaningful, cross-family comparisons less so.

## 4. Structural search run (mind-1's optimizer, my runtime)

`agent/optimize.py` (room event 52 version) on the woven baseline, 12 proposals, six donor
designs, seed 1 (`results/composite/composition_search_fable.jsonl`): 7 of 12 proposals were
geometrically invalid (yinyang donors placed at their library offsets collide with the torus;
duplicates and scale > 0.97 leave the ball), the one valid structural move (replace by a
solenoid) scored 4.1, and the "improvement" that was accepted (scale 0.971, score 20.2 vs 17.5)
is within 160-marker noise. No valid heterogeneous candidate was generated — the same coverage
failure mind-1 saw locally (event 111). Suggested repairs: validate geometry before spending an
evaluation and resample the placement (shrink-to-fit the donor until clearance ≥ 2a), draw
placements from free space around the host rather than from donor offsets, and seed the donor
library with the physics templates that actually fit a torus (rings in the hole, plugs, external
windings) rather than cube-8 offsets.

## 5. Where this leaves the goal

The composition capability is real and usable; the honest confinement result so far is that
**no tested composite beats the bare woven meshed torus** (score 17.2 at 480 markers,
score_plasma equal since it has surfaces and shorting transform), and several of the apparent
gains in the 160-marker screens are evaluator artifacts (A, B). The next things worth an agent's
budget, on evaluator v2.3 with a fixed horizon: (i) l = 2 windings at 3–10× the strand current
(only then is the helical field comparable to the toroidal one; buildability will pay for it),
(ii) placement-aware structural search seeded with torus-compatible templates, (iii) a decision
on the normalisation ROI.
