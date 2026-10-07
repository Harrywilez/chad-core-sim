# The design loop: components → assemblies → score → agent (sessions 7–8, 2026-09-03/05)

*Yes, the system now exists: any design (a list of components, each with a family, its
parameters, and a scale / reflection / rotation / translation / current) can be evaluated
in ~10–30 s by a grading function built strictly from simulation values, and an agent can
iterate on it. Two Opus agents have already run. This file says what was built, what the
score is, what the agents found — including four, then four more, real flaws in the score
that they caught and I fixed — and where it stands.*

## 1. The pieces

| piece | file | what it does |
|---|---|---|
| design spec | `ccsim/design.py` | JSON → `MultiWinding`. 15 component families (`FAMILIES`), transform order scale → reflect → rotate → translate, per-component relative current; multi-circuit families keep their internal senses; closed knots stay closed, open strands get a remote return. `validate()` rejects conductors that touch (clearance < 2 × wire radius) or leave the ball. |
| grading function | `ccsim/evaluate.py` | design → grid (numba Biot–Savart, 2 s) → batched field-line topology → full-orbit Boris (numba, ~1 s per 160 markers × 20 transits) → τ_all, τ_pass, buildability → **score**; optional `--surfaces` → Poincaré / ι → **score_plasma**. Every evaluation is appended to `results/design_ledger.jsonl`. CLI: `python3 -m ccsim.evaluate designs/x.json [--surfaces]`. |
| fast kernels | `ccsim/fastfield.py`, `ccsim/fastorbits.py` | numba Biot–Savart (5e-15 agreement with the reference, 25× faster) and Boris (identical S4 to `run_orbits`, 19× faster); both checked in `exp01_verify`. numpy fallbacks when numba is absent. |
| agent protocol | `agent/PROTOCOL.md` | the objective, the format, the commands, the rules (hypothesis → result per evaluation, never edit `ccsim/`, report flaws rather than exploit them). |
| baseline optimiser | `agent/optimize.py` | random-mutation hill climb for polishing a basin. |
| seeds | `designs/*.json` | the catalogue's best (precess pairs, meshed torus II, mirror pair, solenoid …) and everything the agents made (`agent_*`, `agent2_*`). |
| records | `agent/notes.md`, `agent/leaderboard.md`, `agent/REPORT_run1.md`, `agent/leaderboard_run1_scorev1.md` | hypothesis logs, rankings, the first agent's report verbatim. |

## 2. The score (v2.1; v2.2 changes one clause of score_plasma — §7)

`score = τ_c × buildability`, with `τ_c = √(τ_all · τ_pass)`:

* **τ_all** — confinement time, in wall transits, of 15 keV D⁺ markers seeded uniformly in the
  vessel (r ∈ [0.18, 0.58], ≥ 0.03 from any conductor), isotropic pitch, full orbits on the
  exact vacuum field with the current set so that **B_rms = 2 T** in the region of interest
  (every design at the same field strength). τ = ∫S dt + a fitted exponential tail capped at
  one more run-length; 40 transits, extended to 80 and 120 while ≥ 10 % is still confined.
* **τ_pass** — the same for the markers mirror trapping cannot hold: those inside the loss
  cone of their own field line (ξ₀² ≥ 1 − 1/R_line) plus every marker on a closed line. A
  plasma needs these held (collisions refill the loss cone), so a trap that keeps only
  deeply trapped pitch angles is not rewarded as if it held everything.
* **buildability** = min(1, J_REBCO/J), J = I_max/(πa²), a = the conductor radius the field
  model actually uses (min(design wire radius, 0.45 × clearance)), I_max including the
  largest relative component current.
* **score_plasma** (with `--surfaces`): × 1 if ≥ 20 % of seeds lie on nested surfaces with
  |ι_edge| ≥ max(0.02, (v_d,i/v_th,e)·πR/a) — the polarization is shorted by the electrons —
  else, for closed-line designs, capped at the polarization time τ_pol; open-line designs
  carry no factor (their passing loss is already in τ_pass).

Everything is a simulation value of the design being scored; nothing is fitted to the
catalogue. Compare designs on `score`; rank finalists on `score_plasma`.

## 3. What happened when agents ran it

**Run 1 (Opus, 35 evaluations, score v1).** The agent found that the fat single torus-knot
strand (`hopftorus` η 0.88, 32 circuits, 2 revolutions) scored 309 — twelve times the
meshed torus — and then wrote, correctly, that it did not believe it: every top design had
zero closed field lines and mirror ratio ≈ 1.2, and `--surfaces` found no surfaces. It was
gaming the metric: the fat knot is a dipole-like trap whose *deeply trapped* markers sit on
closed drift orbits forever in vacuum (21 % of the seed population), and v1's tail
extrapolation turned "21 % forever" into τ = 309. Its report listed four flaws — an
inconsistent wire radius between the two buildability terms, the tail lottery, relative
currents ignored in J, and "the score cannot see whether the field confines a plasma" — all
four confirmed and fixed (τ_pass, the tail cap, one conductor radius, I_max with relative
currents). Under v2 the knot scores 6; the closed-line designs lead.

Its physics findings that survived the re-scoring: opening the weave (`delta_eta` 0.03 →
0.08) raises clearance 0.009 → 0.015 and buildability 0.26 → 0.93 at no cost to the field
(`agent_m24w`); splitting a design into several cores loses every time because the seeded
population is what's inside the tube; nested shells are geometrically blocked; the precess
family cannot be made buildable at 2 T.

**Run 2 (Opus, 17 evaluations, score v2).** Pushed the torus fatter (η 0.90: score 20 vs
15.7) and found that at that fatness chirality is nearly free but nearly useless (ι 0.004
at ε 0.15), that an internal ring on the axis intercepts the lines it twists, that shrinking
the torus to fit the wall raises the current it needs, and that a single-turn levitated
dipole needs 2.2 MA for 2 T. It then listed four more flaws: open-line designs escape the
polarization factor; the shorting criterion used the ion thermal speed and was unreachable
(ι ≈ 0.25–0.6 needed) so `score_plasma` had collapsed into a plasma-*size* metric; the tail
cap quantised τ into one bin; the surface radius is seed-range-limited. The second and
third were real and are fixed (electrons carry the shorting current — ι ≳ 0.02 suffices;
the run ladder is now 40/80/120 with the extension keyed on what is still confined); the
first is by design (mirrors are penalised through τ_pass instead) and the fourth is a
measurement limit, both now documented.

**Where it stood at the end of session 7 (score v2.1, geometry later found faulty — see §7–9 for the current state).** The open-weave chiral meshed tori led on both numbers:
`agent2_m24.json` / `agent_m24w.json` (hopfmesh η 0.60, ε 0.70, n 4, δη 0.08, ~5 mm wire)
score **40 / 38**, τ_all ≈ 37 and τ_pass ≈ 50 transits, 47–51 % closed lines, and — with
`--surfaces` — 31 % nested surfaces and ι 0.10–0.16, so they are exempt from the
polarization cap (`score_plasma` 38). The fat untwisted pair (`agent2_r90N24`) scores 18
but is capped to 6 for having no transform; the fat knots score 18–19 with τ_all ≈ 60–70
and τ_pass ≈ 6 — long-lived trapped particles, nothing for the passing ones. The metric now
ranks the design the rungs say would hold a plasma first. Full ranking (every design ever
made, re-scored): `agent/leaderboard.md`.

## 4. What this taught about running agents on a physics score

* The agents were most valuable as **adversarial testers of the score**. Both found the
  exploitable term within a few evaluations, exploited it, and then said so — because the
  protocol asks for a hypothesis per evaluation and a flaw list in the report. Eight flaws
  were found in two runs; six were real. Read the flaw section of an agent report before
  its leaderboard.
* A score built from single-particle vacuum physics is gameable by anything that holds a
  *sub-population*: deeply trapped particles, particles near strands. The fix each time
  was to ask what a plasma needs (passing particles held; polarization shorted), not to add
  an ad hoc penalty.
* Statistical noise matters: 160 markers give ±0.04 on S20, so ordering within ~20 % is not
  resolved; both agents flagged this and asked for 600–1000 markers on finalists.
* The rungs and the score have to agree. Rung 3's original shorting criterion (ions) was
  wrong by √(m_i/m_e); it only surfaced because the agent noticed `surfaces_ok` was never
  true. RUNGS §3 has been corrected.

## 5. How to run the loop yourself

```bash
cd "Chad Core Sim"
python3 -m ccsim.evaluate designs/agent_m24w.json --surfaces   # score + score_plasma of the current leader
python3 agent/optimize.py designs/agent_m24w.json --iters 20   # polish it
# or hand agent/PROTOCOL.md and a budget to an agent (Claude Code / Cowork subagent):
#   "read agent/PROTOCOL.md, then optimise; 40 evaluations; report top 5, lessons, flaws"
```

Next things worth exposing to the agent: the winding-surface deformation (`torus_deformation`
— helical axis, ellipticity, triangularity — is in `geometry` but not in the design JSON),
external helical shaping coils as a component, and a 600-marker re-rank of the top ten.

## 6. Session 8 — watching the designs (the Lab's Designs tab, live)

The Chad Core Lab has a **Designs** tab: every evaluation in `results/design_ledger.jsonl` as a
row (score, τ_all, τ_pass, buildability, S20; hover for closed-line fraction, current, S_end,
time), a score-vs-evaluation chart with the best-so-far line, sort/filter, "one row per design",
and click-to-rebuild: the design JSON is re-generated in the browser with the same 16 component
families and the same scale → reflect → rotate → translate transform as `ccsim.design`
(cross-checked family by family: identical point counts, positions agree to ≤ 1e-12 — the
Hopf-drift family differs by 2e-5 from its arclength table). A "Design from JSON" family in the
geometry panel lets you paste or edit a design and build it; the Readout tab shows the ledger's
score card next to the browser's own field, lines, particles and Poincaré section.

Live, two ways:

* **Locally** — `python3 ui/serve.py` serves the project and the Lab at
  `http://localhost:8765/ui/chad_core_lab_standalone.html`; the Designs tab polls
  `/api/ledger?after=N` every 3 s, so every `ccsim.evaluate` call (yours, `agent/optimize.py`'s,
  or an agent's) appears as it lands, new rows flash, and with "follow the newest evaluation"
  ticked the 3-D view rebuilds to each new design. `/api/status` shows `results/agent_status.json`
  if a run maintains one.
* **In the published Lab** — the same tab subscribes to a shared feed (the artifact's database,
  collection `designs`, one document per evaluation in the compact row format of
  `ui/ledger_rows.py`, plus a `status/run` document). A Claude session pushes rows into it while
  an agent runs (`ui/feed_export.py` writes them as JSON files for the push); the built-in
  snapshot is what you see when nothing has been pushed.

## 7. Session 8 — the geometry was wrong, and fixing it moved the leaderboard

Cross-checking the browser generator against Python turned up two defects in the meshed torus
itself, not in the score:

1. **The weave let the strands touch.** The over/under offset ±δη was set at each crossing and
   ramped between crossings. Where the helical modulation makes a turn of A run nearly parallel to
   a turn of B — the "reversal" bands at ε ≳ 0.5 — the two strands cross twice within a hair, or
   come within a hair without crossing; a ramp through zero there puts both conductors on the
   surface and they meet. At η 0.60, ε 0.70, n 4 the true minimum distance between the strands
   was **0.0013–0.002** (finer polylines find it smaller still), i.e. the conductors overlapped.
   The old clearance estimate sampled every third point across strands and reported 0.015, so the
   score's buildability term never saw it; the agents' best discovery ("open the weave to
   δη 0.08, clearance 0.009 → 0.015, buildability 0.26 → 0.93") was an exploit of that estimator.
2. **Two smaller ones**: the crossing scan could miss the crossing that sits exactly on the
   strand's start (an odd count = one weave defect at the wrap), and the turn-position fixed
   point was iterated 12 times (a 1 % residual at ε 0.7).

The weave is now a **zone weave** (`ccsim.geometry._zone_weave`, ported line for line to the
Lab): wherever the in-surface distance between the strands is below the normal separation
h = 2 δη |∂p/∂η| the two over/under profiles are frozen at opposite full values, ramps happen
only in the free stretches between such zones, and zones facing each other are joined so every
contact has A over and B under; the plain-weave checkerboard sign is kept wherever the strands
are far enough apart to allow it. Clearance is now exact (KD-tree over every point). The
consequence is a regime map that was invisible before:

| δη (η 0.60, ε 0.70, n 4, 24 circuits) | weave | free to weave | min A–B | same-strand spacing | true clearance |
|---|---|---|---|---|---|
| 0.010 | woven | ≈ 50 % | 0.0067 | 0.020 | 0.0067 |
| 0.015 | mixed | 32 % | 0.0100 | 0.020 | 0.0100 |
| 0.030 | nested | 0 % | 0.020 | 0.020 | 0.020 |
| 0.080 | nested | 0 % | 0.049 | 0.020 | 0.020 |

At δη ≳ 0.03 the "meshed torus" is not meshed: strand A sits on the outer torus and B on the
inner one and they never interlock. That changes the *field*: the toroidal currents of the two
strands no longer cancel locally, so the gap between the shells carries the coaxial poloidal field
of the inner strand's net toroidal current (±25 I for 24 circuits), and the interior field is no
longer the clean toroidal solenoid the chiral modulation was designed to twist. Every "open
weave" design lost its closed lines (frac_closed 0.47 → 0.00) and its score (agent2_m24: 40 → 5.5;
agent2_m24n3: 39 → 7). The meshed torus you would build is the **δη ≈ 0.01 woven** one, with a
conductor radius of ~3 mm at 1 m — §8 has what that scores.

The fat untwisted pairs (`agent2_r90N24`…`r100N24`, η 0.90–1.00, δη 0.06 → nested) went the
other way, 17 → 46–57, on the strength of a *slow leak*: not a single closed field line, mirror
ratio 1.1, every line reaches the wall both ways — but the leak is slow, so 160 markers keep 44 %
after 20 transits and the passing markers are "held" for 54 transits. A plasma in such a field
polarizes and leaves on τ_pol (≈ 9 transits here) exactly as in any toroidal field without
transform, and the v2.1 plasma factor missed it because it only capped designs with closed lines.
**Score v2.2** caps every toroidal-field design without nested surfaces and shorting transform —
closed-line or leaky (mirror ratio < 1.5); mirror designs (R ≥ 1.5) still carry no factor. That
is the third exploit found in this project and the first found without an agent: the geometry
change moved the designs into a corner of the score nobody had visited.

`agent/rescore.py` re-evaluates every design (40 min); `agent/leaderboard_from_ledger.py` rebuilds
the ranking from the ledger's latest rows without recomputing. Rows before 2026-09-05 describe the
old geometry.

## 8. Session 8 — "a little suction above and below the hole": tornado plugs

The question was whether a tornado-shaped component above and below the hole could pull escaping
particles back in. Two physics facts first, then the simulation.

* A static magnetic field cannot pull. The Lorentz force is perpendicular to the velocity and does
  no work, and the only average force on a gyrating particle is −μ∇B, which points toward *weak*
  field: a converging field is a *plug*, not a drain — it reflects what approaches it (the mirror
  effect the whole PLAYBOOK is built on). ∇·B = 0 forbids a field that only converges; that would
  be a magnetic monopole. A tornado-shaped field (helical and converging) is perfectly classical —
  it is what a conical helix coil makes, and it is the field of magnetic nozzles and mirror throats
  — what is not classical is net inflow.
* What can pull is an electric field (an electrostatic well: IEC / Polywell virtual cathode) or a
  time-varying magnetic field — the literal rotating tornado, a rotating magnetic field, which
  drives azimuthal current in FRC experiments. `ccsim.fields` has the electrostatic and induced-E
  solvers; neither is in the design loop yet.

So the honest test is the plug. A `cone` family (conical helix, base radius 0.30 → tip 0.06,
height 0.52, 8 turns, tip pointing into the hole from z = ±0.10, base at z = ±0.62) was added to
`ccsim.design` and to the Lab, and the current leader was evaluated with two plugs in the same
sense (a vertical field through the hole), in opposite senses (a cusp at the midplane), with plain
ring plugs, and at relative currents 0.1–1.0. On the old geometry the same-sense tornado plugs at
0.3 scored 44.9 against the bare torus's 40.2 (within the ±10–20 % noise of 160 markers), the cusp
plugs destroyed the closed lines (score 7.8 at I = 1), ring plugs scored 27. On the corrected
geometry every version sits at 4.6–6.4 with the bare torus at 5.5 — the plugs neither help nor
hurt a design whose confinement has already collapsed by nesting. The right host for the plug
test is the woven torus; §9 has that comparison.

## 9. Where it stands after the correction (score v2.2, corrected geometry)

The meshed torus II rebuilt in the woven regime (`designs/woven_*.json`, conductor radius fitted
to 45 % of the exact clearance):

| design | δη | clearance | conductor | τ_all | τ_pass | closed | surfaces · ι (axis → edge) | build | **score** = score_plasma |
|---|---|---|---|---|---|---|---|---|---|
| η 0.60, ε 0.70, n 4, 24 circuits | 0.015 | 0.0100 | 4.5 mm | 36 | 47 | 58 % | 44 % · −0.10 → −0.20 | 0.42 | **17.4** |
| same | 0.010 | 0.0067 | 3.0 mm | 27 | 34 | 59 % | 50 % · −0.10 → −0.25 | 0.19 | 5.9 |
| η 0.50, ε 0.55, n 3 | 0.010 | 0.0082 | 3.7 mm | 14 | 25 | 16 % | 50 % · −0.35 → … | 0.22 | 4.0 |
| 16 circuits, η 0.60, ε 0.70, n 4 | 0.010 | 0.0071 | 3.2 mm | 21 | 25 | 58 % | 38 % · −0.11 → −0.18 | 0.14 | 3.1 |
| pure pair (ε 0), η 0.70 | 0.010 | 0.0057 | 2.6 mm | 19 | 19 | 98 % | 62 % · ι = 0 → capped | 0.16 | 2.9 → 2.2 |

So the physics is intact — τ_c 41 transits, half the seeds on nested surfaces, transform 0.1–0.25,
the pure pair correctly capped for having none — and the number that fell is buildability: a
2 T meshed torus wants ~150 kA in a conductor whose radius the weave limits to 3–4.5 mm, i.e.
2–5 × the REBCO engineering current density. The honest score of the design you would build is
**17 (δη 0.015)**, not 40. The fat nested pairs that top the raw score (57) are capped to 5.6 by
v2.2 (6 % surfaces, ι 0.10 on a single seed, τ_pol 5.6 transits). Full ranking:
`agent/leaderboard.md` (rebuilt from the ledger's post-correction rows).

Tornado plugs on the woven δη 0.015 torus (480 markers, so the ±10 % noise of the 160-marker
screen is gone): bare torus τ_c **40.8** (score 17.2); two same-sense tornado plugs at 0.3 of the
strand current τ_c 18.6 (score 9.1); at 1.0 × τ_c 19.0 (score 9.4); the cusp arrangement 18.1 at
160 markers. The plugs' field threads the hole vertically, adds a non-toroidal component to the
tube field and opens the closed lines (58 % → 36–50 %); the particles follow the opened lines out.
The picture the thunderstorm suggested — funnel the escapers back — cannot be had from a static
coil at all (§8), and the best static approximation to it, a converging plug, costs a factor of two
in confinement on this torus. What would act on escaping particles selectively is an electrostatic
or rotating-field component; those are the next things to add to the component set.

What the agents should do next, on the honest geometry: trade turns for conductor (fewer, fatter
circuits raise the same-strand spacing ceiling of 0.02 and the clearance with it, at the cost of
field ripple), trade tube fatness (η) against the inboard compression that sets the clearance,
and look for a modulation that gives transform with less turn reversal (the reversal bands are
where the strands run parallel and the weave has to float). A run under `agent/PROTOCOL.md` with
40 evaluations would settle the first two in an afternoon.
