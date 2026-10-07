# Fable design search, 2026-09-07 → 09-09: what confines best, honestly

*fable, for John and the room. Goal (GOAL.md): let agents combine different components and test many systems for
confinement. Method: 413 evaluations of 365 distinct assemblies by six Fable research agents working in parallel on
one shared job queue (a cloud workspace plus the Cowork VM on John's Mac), evaluator v2.3 (room event 141) at a
FIXED 120-transit horizon, 160 markers with the matched seeds for screening, 480 markers with held-out seeds
777/4242/999 (and grid 41) for every finalist, `--surfaces` where the plasma factor is meaningful. Every number below
is a simulation value of the design named; the ledger of all 413 evaluations is `results/fable2/rows_all_evaluations.jsonl`,
the full result records of the finalists are in `results/fable2/finalist_results/`, the finalist designs in
`designs/fable2/`, every agent's designs, generators, per-batch hypothesis → result notes and reports are in
`designs/fable-*/`, `notes/fable-*.md`, `reports/fable-*.md`, and the queue tooling in `tools/` (README in
`FABLE_BRIEF.md`).*

## 1. Headline

**A levitated 3-turn ring coil inside an opposing "belt" field is, by the room's own score, an order of magnitude
better than anything the room had — 222 against the woven meshed torus's 15–20 — and it is simple: five circular
coils.** `designs/fable2/magnetosphere_R40tri_belt_i92.json`: three rings of radius 0.40/0.40/0.44 (a compact
triangular bundle, pitch 0.041 ≈ the evaluator's own conductor cap, so the pack is one 1.85 cm-radius conductor × 3
at 1 m scale) carrying +1 each, plus two belt rings of radius 0.90 at z = ±0.25 (outside the wall at r = 0.93) carrying
−0.92 each. At the 2 T-rms normalisation the ring turns carry 835 kA each (J = 0.79 J_REBCO, buildability 1.00).
Every field line inside a spherical separatrix of radius ≈ 0.71 closes around the coil; 92 % of all seeded 15 keV
deuterons — trapped and passing alike — are still confined after 120 wall transits with held-out seeds; τ_all = τ_pass
= 222 transits against the ceiling of 240. It is the magnetic field of a planet: a dipole in an antiparallel external
field ("dipole + uniform field" has an exactly spherical separatrix), i.e. a levitated-dipole experiment (LDX) with its
shaping field chosen to put the separatrix half-way between the seed shell and the wall.

Second, with real rotational transform: the **Levitron** — the fat nested Hopf pair (η 0.85–0.90, 24 circuits) with an
internal ring on the magnetic axis at 2–3 × the strand current — is a plateau at **score 92 ± 4**, buildability 1.0,
true ι from −3 near the ring to −0.2 at the edge of its nest (a tokamak whose plasma current is carried by a levitated
conductor). Third, the plasma-legitimate stellarators without any conductor in the plasma: the **woven chiral torus at
48 circuits** (41.4 held-out, buildability 0.93, 31 % nested surfaces, ι −0.13 → −0.18, score_plasma = score) and a
**helical-axis stellarator made of 30 identical planar circular coils** (34–41 held-out, buildability 1.0, 100 %
nested surfaces, ι 0.31 → 0.42) — both about double the reference woven torus (14.6–20).

## 2. Leaderboard (fixed 120 transits; "screen" = 160 markers matched seeds; "held-out" = 480 markers, seeds 777/4242/999)

| # | family | design (`designs/fable2/…`) | screen | held-out | τ_all / τ_pass | build | closed lines | S120 | plasma factor status |
|---|---|---|---|---|---|---|---|---|---|
| 1 | magnetosphere | `magnetosphere_R40tri_belt_i92` — 3-turn triangle bundle R 0.40 + belt pair (0.90, ±0.25) × −0.92 | 216.7 | **222.4** (grid 41) | 222 / 222 | 1.00 | 1.00 | 0.92 | poloidal-only field: cap not applicable (§5.2); evaluator's `--surfaces` gives 14 or 208 by noise |
| 2 | magnetosphere | `magnetosphere_R40tri_belt_i80` — same, belt × −0.80 | 214.2 | 220.1 | 220 / 220 | 1.00 | 1.00 | 0.91 | same |
| 3 | magnetosphere | `magnetosphere_R45n2_belt_i60` — 2-turn column R 0.45 + belt × −0.60 (four circles) | 210.3 | 216.3 | 216 / 216 | 1.00 | 1.00 | 0.89 | same |
| 4 | magnetosphere | `magnetosphere_R40n2_helmholtz_i50` — 2-turn R 0.40 + Helmholtz pair (0.80, ±0.40) × −0.50 | 178.5 | — | 207 / 207 | 0.86 | 1.00 | 0.84 | phase-1 leader; oblate separatrix (§4.1) |
| 5 | dipole | `dipole_ring_R40_helmholtz_i25` — lone ring + Helmholtz × −0.25 | 96.8 | 102.7 (103.6 grid 41) | 211 / 211 | 0.49 | 1.00 | 0.87 | lone-ring conductor cap |
| 6 | dipole | `dipole_R40_sphere_m003` — lone ring + `sphere` winding × −0.03 | 85.7 | 90.9 | 182 / 193 | 0.49 | 0.66 | 0.73 | same physics, winding inside the wall |
| 7 | Levitron | `levitron_e090_N24_R40_I2` — nested pair η 0.90 N 24 δη 0.06 + axis ring R 0.40 × 2.0 | 89.9 | 94.6 (92.9 second run) | 92 / 97 | 1.00 | 0.41 | 0.35 | true ι −1.6 → −0.21; evaluator's ι aliased → capped to 15.6 (§5.1) |
| 8 | Levitron | `levitron_ringdominated_dip40_pair80_040` — ring R 0.40 × 1 + nested pair η 0.80 × 0.4 | 87.9 | 94.3 | 95 / 94 | 1.00 | 0.34 | 0.35 | true ι −4 → −0.3; aliased → 16 |
| 9 | Levitron | `levitron_e085_N24_R40_I3` — η 0.85, ring × 3.0 | 94.7 | 90.3 | 88 / 92 | 1.00 | 0.47 | 0.33 | same |
| 10 | dipole | `dipole_stack_R35_n3` — 3-turn stack R 0.35, no shaping field | 78.7 | 76.9 | 89 / 89 | 0.87 | 0.52 | 0.36 | poloidal-only |
| 11 | fat pair (control) | `fatpair_e100_capwire_control` — nested pair η 1.0, wire at cap | 58.2 | 49.6 | 57 / 43 | 1.00 | 0.00 | 0.17 | **capped ≈ 6**: leaky toroidal solenoid, no transform |
| 12 | woven stellarator | `woven_N48_de016` — hopfmesh η 0.60 ε 0.70 n 4, 48 circuits, δη 0.016 | 45.9 | **41.4** | 38 / 52 | 0.93 | 0.35 | 0.13 | **score_plasma 41.4**: 31 % surfaces, ι −0.13 → −0.18 |
| 13 | woven stellarator | `woven_N24_de020` — 24 circuits, δη 0.020 | 38.9 | 41.3 (screen with surfaces) | 47 / 64 | 0.75 | 0.53 | 0.16 | score_plasma 41.3: 44 % surfaces, ι −0.10 → −0.22 |
| 14 | woven stellarator | `woven_N24_de025_e65_vfp30` — η 0.65 + vertical-field ring pair | 59.0 | 39.1 | 35 / 44 | 1.00 | 0.25 | 0.10 | score_plasma 39.1: 31 % surfaces |
| 15 | planar-coil stellarator | `heliax_Rc52_N30_rk24_M3_d12` — 30 planar circles r 0.24 on a 3-period helical axis | 48.6 | 41.2 (matched 480) | 34 / 50 | 1.00 | 0.41 | 0.11 | score_plasma = score (100 % surfaces, ι 0.36 → 0.42 on the baseline) |
| 16 | planar-coil stellarator | `heliax_Rc54_N30_rk24_M3_d12` | 54.9 | 34.0 | 26 / 44 | 1.00 | 0.26 | 0.08 | score_plasma 54.1 at 160 markers |
| 17 | CNT-type | `cnt_a50_d10_t85_pf95z30_im70` — two interlocked circles + PF pair | 23.9 | — | 35 / 16 | 1.00 | 0.12 | 0.07 | best of its family; cannot compete |
| ref | woven torus | `designs/composite/woven_m24_de015.json` (24 circuits, δη 0.015) | 20.0 | 14.6 | 31 / 39 | 0.42 | 0.59 | 0.09 | score_plasma = score |

The 160-marker screen is optimistic for every toroidal family (the woven and helical-axis designs read 10–40 % high;
the magnetosphere designs, whose losses are ~20 markers out of 160, confirm *upward*). Rank by the held-out column.

## 3. What was searched (six directions, two phases)

Each agent worked from a stated hypothesis per batch, at the fixed horizon, with `wire_radius` set to 0.45 × the exact
clearance (so buildability is judged at the conductor the geometry allows), and with the evaluator's known artifacts
(ROI-normalisation leak, run ladder, lone-conductor chord clearance) checked rather than exploited. Per-agent reports:

* **fable-rings** (axisymmetric internal conductors → the magnetosphere; 92 evaluations, `reports/fable-rings.md`).
  Lone ring R 0.25–0.50; multi-turn coil packs as `circle` stacks (the `solenoid` family is useless here: its remote
  return leads run through the ball and break axisymmetry); same-sense ring pairs; octupoles (bad); ring + opposing
  Helmholtz pair (H7); then in phase 2 the *shape* of the separatrix: Helmholtz → belt, bundle geometry, current ratio,
  ring z-offset, polar pairs, loss anatomy (§4.1).
* **fable-levitron** (fat toroidal-field pair + internal axis ring; 63 evaluations, `reports/fable-levitron.md`).
  Ring current 0.3–9, ring radius 0.36–0.56, η 0.80–1.00, circuits 16–32, δη 0.04–0.10, two rings, ring z-offset,
  scaled tube, opposing/same-sense vertical field pairs, unaliased transform with `tools/iota.py`.
* **fable-weave** (the woven chiral meshed torus made buildable; 80 evaluations, `reports/fable-weave.md`).
  Circuits 10–48 × δη 0.012–0.05, η 0.60–0.80, ε 0.55–0.85, periods 3–6, helical mode, vertical-field ring pairs;
  the weave/nest boundary, and which combinations keep ≥ 20 % nested surfaces.
* **fable-cnt** (stellarator surfaces from planar circular coils; 49 evaluations, `reports/fable-cnt.md`).
  CNT-type interlocked coils (tilt, offset, PF coils both signs), heliacs with a central ring (negative), and the
  helical-axis set of 24–40 planar coils (period number, axis excursion, coil radius and clearance).
* **fable-hybrid** (compositions and a placement-aware structural search; 99 evaluations, `reports/fable-hybrid.md`).
  Phase 1: fat pair + cones / vertical-field rings / picket; dipole + toroidal field / sphere winding / PF rings;
  link, baseball, yin-yang around tori. Phase 2: polar plugs (cones, rings) on the magnetosphere, a toroidal field on
  the magnetosphere, PF shells around the Levitron, split stacks; and `agent/composition.py`'s `propose()` driven from
  the two leaders with a placed donor library, 144 proposals → 102 valid (71 %; phase 1 of the room: 5/12) → 22
  evaluated → 0 improvements.
* **fable-stellar** (external l = 2 helical windings on nested pairs) never got past its first batch in phase 1; the
  woven agent's findings cover the transform-for-buildability trade instead.

## 4. What was learned

### 4.1 The magnetosphere: why it wins, and where its last 8 % go

A ring current plus an antiparallel uniform field has a spherical separatrix of radius r_s = (μ₀ m / 2π B)^{1/3};
inside it every line is closed around the coil, and in an axisymmetric field the canonical angular momentum p_φ
confines every particle whose drift shell stays off the wall and off the conductor — passing particles included, so
τ_pass = τ_all and the passing-particle term that kills mirrors does nothing here. The seed shell (r ≤ 0.58) fits
inside r_s ≈ 0.71 with the wall at 0.82.

* The phase-1 Helmholtz pair (0.80, ±0.40) made an **oblate** separatrix (a real ring's on-axis field is weaker than a
  point dipole's, its equatorial far field stronger): the equator touched the wall while the polar null sat only 0.05
  above the seed shell. Reconstructing the seeds of every job and computing ψ = ρA_φ, *every* wall-lost marker had
  ψ/ψ_max < 0.022 — the polar column. A **belt** pair (ρ 0.90, z ±0.25: more opposing field at the outboard equator,
  less on the axis) makes the separatrix round (z_null 0.71 = ρ_s 0.72), and all nine belt designs beat all five
  Helmholtz ones by 10–15 points.
* The shaping current has a sharp threshold on equatorial closure and then a plateau; the optimum for the R 0.40
  triangle is 0.90–0.92 × the turn current; ring z-offset is strongly bad (+0.05 → 166); coil shape (R 0.35–0.45,
  column / triangle / 2×2 / diamond, 2–4 turns) is a 205–217 plateau once the belt is there, the compact triangle best
  by ~5 because fewest field lines thread the gaps between turns. Buildability was never the lever after the first
  multi-turn step (J/J_REBCO 0.6–0.87 at build 1.0).
* Two floors remain, ~6 % each: the **polar column** (the same nine seeds at ρ ≤ 0.15, |z| 0.37–0.52 die in every
  design — there ψ is below the p_φ excursion m v ρ/q, and holding them would need > 0.5 T on the axis in the hole,
  which no ring gives at the 2 T-rms normalisation) and the **conductor** (seeds within 0.04–0.09 of the coil; for
  bundles the channel is real gap-threading, grid-independent; a lone ring's hits halve at grid 41 — interpolation).
  Polar plugs (cones or rings at the poles, either sign, 0.002–0.3 of the ring current) do nothing or open the
  separatrix — the on-axis null is topologically protected for any axisymmetric coil set; a toroidal field on top of it
  is a cliff (τ 207 → 38 at 5 % of the ring current: ripple and strand interception). The ceiling of the family at
  this normalisation is ≈ 225; 222 is where it stands.

### 4.2 The Levitron: transform for real, at a fixed price

The internal ring must carry ~2–3 × the strand current for the ring's nest of surfaces to hold anything; below that the
old "axis ring I = 1" design (29.6) sits in the wrong regime. The nest has radius 0.19–0.22 around the ring (set by the
pair's residual vertical field, not by the wall), holds ≈ 45 % of the seeds, and its true transform is |ι| 1.3–3.3
near the ring falling to 0.2–0.3 at the edge — far above the 0.02 shorting threshold. Losses: 35/160 to the wall
(seeds in the hole and polar caps, prompt) and 75 "wire", of which 54 are the *inner strand shell* (∇B-drift
excursions of nest-edge markers into the top/bottom of the inner torus), not the ring. Hence the plateau: η 0.85–0.90,
R 0.38–0.44, ring 2–3 ×, N 20–28 all give 86–95; δη 0.04/0.05, split rings, z-offsets, scaled tubes and vertical
fields of either sign are all worse — the perfect nest would give ≈ 110. Levitron (transform, 92) and magnetosphere (no
transform, 222) are the two ends of one trade-off: the toroidal field buys shorting and costs confinement.

### 4.3 The woven chiral torus: go to MORE circuits, not fewer

Clearance is the over/under separation ≈ 0.66 δη at every circuit count; the same-strand ceiling ≈ 0.50/N only binds
at N ≤ 14. So buildability ∝ N δη² and the buildable optimum is **N 40–48 at δη 0.016–0.018** (τ_c 42–49 at build
0.93–1.0), not the "fewer, fatter circuits" of AGENT.md §9 — at fixed δη fewer circuits only raise the current per turn.
The surfaces come in two discrete states set by the zone weave (big: axis R 0.485, surfaces ρ 0.33–0.66, when N·δη ≲ 0.5;
small: axis 0.39–0.42), and the 20 % surface threshold is effectively 4 of 16 Poincaré seeds — several good designs sit
one seed from a ×0.1 cliff. Fatness η ≥ 0.65 raises raw τ but pulls the axis inboard and shrinks the surfaces; ε 0.85
doubles ι at the same surface count; a vertical field from a ring pair moves the axis outboard (+0.3 → more surfaces) but
did not survive the held-out check. Honest family number: **41 at buildability 0.93 with score_plasma = score** (the
reference: 14.6–20 at 0.42).

### 4.4 A stellarator from planar circles

Thirty identical planar circular coils (r 0.24) centred on the helix ρ = 0.54 + 0.12 cos 3φ, z = 0.12 sin 3φ with
their normals along it make nested surfaces (100 % of the Poincaré seeds, ι 0.31–0.42, no islands) with no conductor in
the plasma and no weave; buildability 1.00 at the chord cap. It is limited by coil–coil clearance (r 0.26 gives more
surfaces but build 0.2–0.8) and by the wall cutting the outboard tube (R_c 0.56 collapses). 34–41 held-out. The
CNT-type interlocked pair needs a *negative* PF current and a small offset and confines ≈ 10 % of the shell (23.9);
a heliac with its central ring is worse than the same coils without it (the ring destroys their surfaces).

### 4.5 Compositions: the search itself

Two compositions made the whole difference — ring + shaping field, and torus + axis ring — and both are two-family
assemblies found by physics, not by mutation. Everything added to a converged leader afterwards was neutral or worse
(cones, polar rings, PF shells, picket cusps, a toroidal field on the magnetosphere, split stacks, third belt rings).
The structural search now produces valid heterogeneous candidates (71 % valid with placed donors and a
validate-before-evaluate filter, against 5/12 in the room's phase 1) but none improved a leader: the leaders sit on
symmetric ridges and `composition.py` has no symmetric-pair operator, and a coarse-grid closed-fraction pre-screen ranks
candidates well enough to spend evaluations only on the plausible ones.

## 5. Evaluator findings (documented, not exploited; nothing in `ccsim/` was changed)

1. **Rotation-number aliasing.** `ccsim.surfaces.classify` measures ι from the poloidal angle between successive Poincaré
   punctures with `np.unwrap`, so |ι| > 0.5 is reported modulo 1; the "edge" value is also taken from the outermost
   surviving seed, which for the Levitron lies on near-wall solenoid lines outside the nest. Result: designs with true
   ι −3 → −0.2 get ι_edge ≈ 0.001, "not shorted", score_plasma ≈ 16 at score ≈ 90. `tools/iota.py` integrates θ and φ
   continuously and gives the true profile (30 s).
2. **Purely poloidal fields** never advance in φ; the Poincaré classifier then reports either "surfaces" with a spurious
   ι whose sign is grid noise (208 = uncapped) or none (14 = capped) for the same topology. Physically an axisymmetric
   dipole plasma does not polarize — the ∇B and curvature drifts are toroidal, in the symmetry direction, so no charge
   separation builds up (LDX, planetary magnetospheres) — and the rung-3 cap is a model limitation for this family; the
   real caveats are different: a levitated superconducting coil inside the plasma (levitation, cryostat, particle sink —
   our 6–8 % conductor channel), the interchange pressure-profile limit p ∝ V^{−5/3}, the polar X-points feeding the wall,
   and the vacuum normalisation ignoring high-β diamagnetism.
3. **The surface threshold is 4 of 16 seeds**, a ×0.1 cliff one seed wide; the axis finder's fan starts at R0 − 0.55 r_tube,
   marginal for fat tubes.
4. **The seed filter over-estimates distance near a wire** (trilinear interpolation of a cone), so 2–4 % of markers start
   at a conductor; the 25³ grid (dx 0.07) cannot resolve a ring's 1/d field within ~0.1 (a lone ring's conductor hits halve
   at grid 41; bundles' do not).
5. **τ_pass bookkeeping**: a component that mirror-traps markers out of the passing set (polar cones) raises the score
   with no change in retention — compare retained counts, not the labels.
6. **`solenoid` is not a coil here**: its remote return runs radially through the ball at one azimuth and breaks
   axisymmetry (closed 0.14–0.29 at build 1.0). Multi-turn coils must be `circle` stacks.
7. **Global `wire_radius`**: one conductor radius for all components ties an auxiliary's clearance to the main coil's
   buildability (a nested `hopfmirror` shell at η 1.0 has clearance 0.016 at any N).
8. Noise, as advertised: 160 markers ±10–20 % (worse for toroidal families whose long-lived population is small),
   cross-machine chaos 1–7 %, grid 25 → 41 ±2 points; the 40/80/120 ladder was avoided throughout by the fixed horizon.
9. `propose()` jitter pushes donors at ρ ≳ 0.9 out of the ball (~15 % of proposals) and `duplicate` flips the current sign
   at random; the earlier "7 of 12 invalid" was mostly that.

## 6. What next

* If the room wants score_plasma to rank these designs it needs (a) an unaliased ι (the `tools/iota.py` method, or the
  puncture count per poloidal turn), (b) the edge taken on the outermost *nested* surface, (c) a "poloidal-only" branch
  that does not apply the polarization cap to axisymmetric fields — or an explicit decision that internal-conductor
  designs are out of scope. All three re-rank finalists but do not touch the raw score.
* The magnetosphere at 222 is at its floor under this normalisation; the honest next questions are engineering ones
  (a real levitated coil pack, its cryostat and support, the polar divertor region) and plasma ones the score cannot see
  (the interchange profile limit, high-β equilibrium). The Levitron's ceiling (~110) is set by the inner strand shell:
  a toroidal-field coil that leaves the polar caps of the tube open (a partial shell) would move it.
* For the meshed torus John likes: the buildable version is 48 circuits at δη 0.016 (conductor 4.7 mm at 1 m), 41 with
  surfaces and transform, twice the reference; the search over ε, periods, η and vertical fields is exhausted at 160-marker
  resolution — the next gain there needs 480-marker screening.
* John's request in event 152 (a viewable layer): every design here is a plain component list the Chad Core Lab can
  rebuild in 3-D from its JSON ("Design from JSON" family, Designs tab); the belt design is five circles. The rows in
  `results/fable2/rows_all_evaluations.jsonl` are in the queue's compact format; `tools/qlib.py` documents it.
