# Design-agent run 1 (Opus, 2026-09-02) — report as returned, verbatim

Score version: v1 (τ_all with uncapped tail; buildability = clearance/3a × J term with inconsistent radii).
The flaws it found (F1–F4) were fixed in score v2 (see PROTOCOL.md); the leaderboard was re-scored.

Budget spent: 35 scored evaluations (~30 min of evaluation wall-clock), plus ~40 cheap pre-screens (`agent/prescreen.py` — geometry, and where useful the Biot–Savart grid only; no particle run, no score) used to reject un-buildable geometry before spending an evaluation.

## Top 5 designs (score v1)

| # | score | τ_c | build | clearance | I(2 T) | J/J_REBCO | S4 / S20 / S120 | file |
|---|-------|-----|-------|-----------|--------|-----------|-----------------|------|
| 1 | 309.6 | 309.7 | 1.00 | 0.0211 | 130 kA | 0.46 | 0.38 / 0.24 / 0.21 | designs/agent_knot_e0.88_N32.json |
| 2 | 300.2 | 300.5 | 1.00 | 0.0177 | 104 kA | 0.52 | 0.39 / 0.23 / 0.21 | designs/agent_knot_e0.88_N40.json |
| 3 | 275.0 | 275.1 | 1.00 | 0.0211 | 115 kA | 0.40 | 0.34 / 0.23 / 0.19 | designs/agent_knot_e0.85_N36.json |
| 4 | 219.6 | 296.8 | 0.74 | 0.0155 | 207 kA | 1.35 | 0.69 / 0.31 / 0.20 | designs/agent_knot_e1.15_N20.json |
| 5 | 188.9 | 189.0 | 1.00 | 0.0219 | 262 kA | 0.86 | 0.62 / 0.36 / 0.21 | designs/agent_knot_e1.12_N16.json |

All five are the same object: one closed (N+1, N−1) torus-knot conductor on a fat Hopf torus (`hopftorus`, `revolutions = 2`) — no mirror partner, no weave. Previous best known: 24.1.

Best design in the mesh (stellarator-weave) family: designs/agent_m24w.json, score 62.7 (τ_c 67.3, build 0.93) — hopfmesh η 0.60, ε 0.70, n 4, N 24, `delta_eta` 0.08, wire 0.0049. 2.6× the previous best, and the one I would actually build (see F4).

## The three most useful things learned

1. Buildability is bought with conductor spacing, and spacing is a design variable. J = I/(π a²) with a = min(1 cm, 0.45·clearance), so J ∝ I/clearance². Opening the mesh weave (`delta_eta` 0.03 → 0.08 at N=24) raised clearance 0.0090 → 0.0148 and cut J 2.7×: buildability 0.264 → 0.93, score 24.1 → 62.7 (agent_m24w), with no change to the confining physics. Since I ∝ 1/N and clearance ∝ 1/N, J ∝ N — fewer, better-separated strands are more buildable, up to clearance ≈ 0.022 where a saturates at 1 cm.
2. Dropping the mirror partner beats any parameter tuning, and fatter wins far past the playbook values. A single strand has 2–4× the spacing of the woven pair: agent_knot32 (η 0.70, N 32) scored 58.7 on its first evaluation. Then η 0.75 → 0.88 at N = 32 took the score 102.7 → 309.6, and the survival curve kept improving to η ≈ 1.15 (S4 0.67, S20 0.43, S120 0.26). Mechanism: retention = fraction of the seeded shell inside the conductor's tube, so anything that grows the tube wins.
3. Splitting into several cores loses every time, for that same reason. Twin z-mirrored knots (0.80 each, z = ±0.40) scored 5.7 (same sense) and 8.5 (opposed). Shrinking one mesh to scale 0.85 did the same (closed fraction 0.51 → 0.27, score 62.7 → 7.0). Nested shells are geometrically blocked (clearance ≤ 0.005 for every η pair tried). Also falsified: same-sense currents (closed fraction 0.00, τ_c 6.4) and a vertical-field circle pair on the mesh (score 27.2, survival unchanged). The precess family was rejected at pre-screen: its best full-scale clearance is 0.0070, so J/J_REBCO ≳ 5.

## What I would try next

- Expose the winding-surface deformation (torus_deformation) from the design JSON.
- Odd revolutions / other knot classes (revolutions 6, 10 change the toroidal/poloidal current ratio — the transform).
- A finer (ε, delta_eta) mesh scan at N = 24.
- Re-rank the top ten with ~1000 markers before believing any ordering inside it.

## Suspected flaws in the grading function (all four confirmed and fixed in score v2)

- F1 — the two buildability terms used different wire radii (declared vs used); declaring a thin wire set one term to 1 for free.
- F2 — τ_c was a lottery: ±2× generally, ±10× at the ladder boundary; the uncapped tail term dominated (knot η0.75 N40 scored 10.1 with S20 = 0.16 while knot η0.75 N32 scored 102.7 with S20 = 0.15).
- F3 — I_max ignored relative component currents in J.
- F4 — the score could not see whether the field confines a plasma: every design in the top 18 had closed-line fraction 0.00, connection length 2–5, mirror ratio 1.1–1.35, and `--surfaces` on the best knot reported 0 % flux surfaces, ι undefined. "I do not believe the knot is a better confinement device than the meshed torus — I believe it is better at the metric."

Recommendation (v1): for the number, designs/agent_knot_e0.88_N32.json; for a plasma, designs/agent_m24w.json — the open-weave chiral meshed torus, the only design above 60 with closed field lines and a family whose surfaces and ι ≈ 0.1–0.24 are verified in RUNGS §1.
