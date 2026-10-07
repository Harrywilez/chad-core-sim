# Would the meshed torus actually work? — the fidelity rungs (session 6, 2026-09-02)

*Rung 1 (flux surfaces), rung 2 (guiding-centre transport), rung 3 (collective
polarization loss, as bookkeeping) and rung 4 (strand forces and conductor limits) for the
woven Hopf mirror pair.  Reproduce with `experiments/exp10_surfaces.py`, `exp11_transport.py`,
`exp12_rung_figures.py`; numbers in `results/exp10_*.json`, `exp11_*.json`,
`exp12_forces.json`; the Lab has the new "Meshed torus II" family, a Poincaré-section button,
and rung-1/rung-2 tables in the Ledger.*

## 0. The answer

> **Session-8 note (2026-09-05).** The rung-1 and rung-2 numbers below were computed with a weave
> whose two strands touched in places (AGENT.md §7) at δη 0.025–0.03. The weave has been rebuilt
> (zone weave, exact clearance) and re-scanned in §7: the physics of this section — chiral
> modulation gives nested surfaces with ι 0.1–0.25, the pure pair has none — holds, *provided the
> torus is actually woven*, i.e. δη ≲ 0.015 (conductor radius ≲ 3–5 mm at 1 m). At δη ≳ 0.03 the
> two strands nest instead of interlocking and the surfaces disappear. Read §7 before quoting a
> number from §0–§2.

**As drawn, no. Modified in one specific way, it becomes a real (weak) stellarator — and the
modification is a change to the weave, not to the shape.**

The single-particle screen (PLAYBOOK.md) made the meshed torus look like the winner because
it confines what starts inside the tube. Rung 1 shows why that cannot survive contact with a
plasma: the pair is **achiral** — reflecting it and rotating it by 180° gives it back with
its currents reversed — and an achiral field has **zero rotational transform**. It is a
toroidal-field coil and nothing else. Every field line is a circle that any imperfection
turns into a slow outward spiral (a 0.03 % stray radial field opens the lines within 30
turns), and, more importantly, a plasma in a pure toroidal field polarizes and E×B-drifts
out as a body on ~50 transits (≈ 30 µs at reactor scale, rung 3). No amount of twisting the
*torus* fixes this: a poloidal current sheet on a rotating-ellipse or helical-axis winding
surface is still a solenoid inside (ι ≈ 0.02 in the scans).

What does fix it is making the **current pattern** chiral. If the two strand families are
woven with *opposite* helical density modulations, 1 ± ε cos(2θ − nζ), the poloidal current
stays uniform (still a toroidal solenoid) while the toroidal current acquires the
cos(2θ − nζ) pattern of a classical l = 2 stellarator's helical windings. That field has
nested flux surfaces with real transform:

| chiral meshed torus (η, ε, n) | surfaces | last surface r | ι axis → edge | well |
|---|---|---|---|---|
| η 0.50, ε 0.55, n 3 | 42 % of seeds | 0.122 | 0.36 → 0.42 | −0.09 |
| η 0.50, ε 0.70, n 4 | 33 % | 0.126 | 0.32 → 0.42 | −0.09 |
| η 0.60, ε 0.85, n 4 | 42 % | 0.138 | 0.15 → 0.33 | −0.21 |
| η 0.60, ε 0.70, n 3 | 33 % | 0.115 | 0.21 → 0.29 | −0.08 |
| η 0.60, ε 0.70, n 4 | 46 % | 0.158 | 0.10 → 0.24 | −0.19 |
| η 0.60, ε 0.70, n 5 | 50 % | 0.171 | 0.05 → 0.18 | −0.24 |
| η 0.70, ε 0.70, n 4 | 42 % | 0.170 | 0.03 → 0.10 | −0.17 |
| η 0.60, ε 0 (pure mirror pair) | 8 % | 0.013 | 0 | — |

(unit ball; the tube radius is 0.36–0.39, so the confined region is 30–45 % of the tube
radius. Poincaré sections: fig12; the scan: fig13.)

Rung 2 then confirms the physics on particles: 15 keV deuterons seeded inside the surface
region of the untwisted torus are 97 % still there after 40 transits and **58 %** after 60 —
the vertical drift walks them out on schedule; in the chiral torus (η 0.60, ε 0.70, n 4) it is
95 % and **88 %** (fig14). Rung 3's shorting criterion (any transform above ~0.02, since the
electrons carry the shorting current) is met by the chiral configurations, not by the pure
pair.

So the object you would build is: two closed (N+1, N−1) torus-knot strands, exact rotated
mirror images, woven, opposed currents, with the weave density bunched into n helical bands
in opposite senses on the two strands. Visually it is still the meshed torus (fig11 vs the
Lab's "Meshed torus II" family); electrically it is a stellarator.

## 1. Rung 1 — flux surfaces

### 1.1 What was built

`ccsim/surfaces.py`: batched RK4 field-line tracing on the **exact** Biot–Savart field
(`ccsim/fastfield.py`, a numba kernel that reproduces `fields.biot_savart` to 5e-15 and is
25× faster), Poincaré punctures of the φ = 0 half-plane, magnetic-axis search (coarse fan →
smallest puncture spread → centroid, refined), rotational transform per line, Fourier
surface test (m ≤ 6, residual < 6 %), specific volume ∮dl/B per turn and the well depth
1 − U_edge/U_axis. Grids are deliberately not used: trilinear interpolation is not
divergence-free and would manufacture or destroy surfaces.

### 1.2 Two defects in the original object, fixed

* **Open strands.** The 24 × 1 strand ends one poloidal turn short of closing, and that
  missing turn leaves a φ-averaged radial field ∝ sin(τ₀)/N at the midplane — 1.2 % of B_tor
  for the session-5 geometry. With 24 circuits × **2 revolutions** each strand is a closed
  (25, 23) torus knot: no leads, no gap, and the pair is exactly inversion-symmetric.
* **The weave itself was chiral by accident.** Displacing A outward and B inward at
  alternate crossings is a checkerboard only if the over/under pattern is consistent with B
  being the exact mirror of A. The old (−1)^k pattern was not (B carried its own profile),
  leaving a ~1 % stray B_z; the fix is to alternate by *rank along the strand* — for
  (N even, 2 revolutions, mirror rotation π) every crossing pair has odd rank sum, so one
  profile serves both strands and B = R_π M A exactly. Stray field at mid-radius drops from
  0.6 % to 0.06 %; the weave is still 99.6 % plain.

### 1.3 The achirality theorem, and what it rules out

With B = R_π·M·A and I_B = −I_A the pair maps to itself under inversion with all currents
reversed, so **B(−x) = −B(x)**. Inversion maps a field line to a field line of the opposite
handedness through the image point; nested surfaces are inversion-symmetric sets; hence
ι = −ι = 0 on every surface. That is exact, and it is why every "shape" deformation failed:

| deformation of the winding surface (η 0.70, opposed, 24 × 2) | surfaces | ι |
|---|---|---|
| none | 38 % (toroidal circles) | 0 |
| rotating ellipse ε 0.10 / 0.20 / 0.30, n 5 | 44 / 44 / 12 % | 0.02 |
| rotating ellipse 0.20 with n 3 / 7 | 38 / 12 % | 0.02 / 0.02 |
| helical axis 0.15 / 0.30 | 19 / 0 % | 0.02 / — |
| ellipse 0.20 + helical axis 0.15 / 0.30 | 0 % | — |

The residual 0.02 is the small symmetry breaking of the deformation itself. A poloidal
current sheet has a uniform toroidal field inside whatever its cross-section, like a
solenoid of any shape — twisting the sheet does not twist the field.

### 1.4 Chirality from the weave

Two ways to modulate the strand density were tried:

* **common** modulation (both families 1 + ε cos): the *poloidal* current carries the
  pattern → a toroidal-field ripple, mirror-like, ι ≤ 0.02.
* **opposite** modulation (A: 1 + ε cos, B: 1 − ε cos): the poloidal current stays uniform
  and the *toroidal* current becomes 2ε cos(2θ − nζ) — the l = 2 helical-winding current.
  This is the one that works (§0 table).

A subtlety that cost an afternoon: modulating the local pitch of a strand does nothing,
because the modulation averages out over each toroidal transit. The turns must be
*displaced*: they are the level sets of a2 + (ε/2) sin(2a2 − n a1) − (q/p) a1, so the turn
density is 1 + ε cos(2a2 − n a1) exactly (`geometry._modulated_strand`). Crossings of the two
modulated families are then found exactly in the (a1, a2) torus coordinates and woven by
rank (`geometry.hopf_helical_pair`, 87 % plain alternation along the second strand).

Trends (fig13): ι ∝ ε² roughly; thinner tube (smaller η) gives much more transform for the
same modulation (η 0.70 → 0.60 → 0.50 at ε 0.70, n 4: 0.10 → 0.24 → 0.42) at the price of a
smaller plasma; fewer periods (n 3–4) beat n 5–7; ε beyond ~0.85 crowds the strands into
each other (clearance < 0.005) and the good-surface region shrinks. Almost every chiral
case sits in a magnetic **hill** (well depth −0.1 to −0.2) — interchange-unstable at finite
β without further shaping — the one thing in this table that a real design would have to
fix next (a helical-axis excursion or triangularity on top of the modulation; the vacuum
machinery to test it is all here).

The magnetic axis also moves a long way inboard with the modulation (R 0.64 → 0.48 for
η 0.60, ε 0.70, n 4): seeding from the geometric centre found nothing, which is why the
analysis now locates the axis first.

## 2. Rung 2 — guiding-centre transport

`ccsim/guiding_center.py`: RK4 on the guiding-centre equations in a vacuum field
(dX/dt = v∥b + (v∥² + v⊥²/2)/Ω · b×∇B/B, dv∥/dt = −μ b·∇B/m, μ conserved), B and ∇B from
the exact kernel by central differences, Lorentz pitch-angle scattering at fixed speed with
an adjustable rate ν. Loss = the sphere r = 1.05 or 0.02 from a strand (the winding is the
vessel). Validation: energy is conserved to round-off with ν = 0, and particles seeded
inside the tube of the untwisted torus are 100 % retained for 10 transits with or without
ν = 2/transit scattering — the toroidal field has no loss cone to scatter into.

| run (15 keV D⁺, 64 particles, seeded inside 0.7 × last-surface radius) | S5 | S20 | S40 | S60 |
|---|---|---|---|---|
| untwisted η 0.60, ν = 0 | 100 % | 100 % | 97 % | **58 %** |
| untwisted η 0.60, ν = 0.1 / transit | 100 % | 100 % | 98 % | 72 % |
| chiral η 0.60, ε 0.70, n 4, ν = 0 | 100 % | 100 % | 95 % | **88 %** |
| chiral η 0.60, ε 0.70, n 4, ν = 0.1 / transit | 100 % | 100 % | 91 % | 84 % |
| chiral η 0.60, ε 0.70, n 4, ν = 0.5 / transit | 100 % | 95 % | 78 % | 53 % |
| chiral η 0.50, ε 0.55, n 3, ν = 0 | 97 % | 92 % | 86 % | 83 % |

The untwisted curve is the drift: nothing happens for 40 transits and then the population
reaches the strands together (fig14), exactly the v_d ≈ 2T/(eBR) ≈ 9 km/s of rung 3
(0.25 m in ~28 µs ≈ 40 transits); scattering does not change that, because in a toroidal
field the drift direction is set by the charge, not the pitch angle. The chiral torus
loses 12 % by 60 transits with ν = 0 — helically trapped particles, the stellarator's
classic collisionless loss channel, which a real design reduces by quasi-symmetry
optimisation — and the loss grows with collisions (16 % at ν = 0.1, 47 % at ν = 0.5 per
transit) as scattering keeps feeding the trapped/loss orbits. The higher-transform,
thinner-tube case (η 0.50, n 3) has a prompt 8 % loss in the first five transits and then a
flatter tail (83 % at 60): more transform, but a stronger helical ripple to be trapped in.
On the scale of ν: the physical pitch-angle rate for 15 keV D⁺ in a 10²⁰ m⁻³, 10 keV plasma
is ~10⁻⁴ per transit, so the ν = 0.1–0.5 runs are 10³–10⁴× that; at the physical rate the
collisionless 0.2 %/transit ripple loss is what matters, and it is the number to optimise.

## 3. Rung 3 — the collective loss, as bookkeeping

`ccsim/polarization.py`. In a pure toroidal field the ion and electron drifts separate
charge vertically; the resulting E×B drift moves the plasma out as a body at ≈ v_d:

| configuration | v_d (15 keV, m/s) | τ_pol (transits) | ι needed to short | ι available |
|---|---|---|---|---|
| untwisted, R 0.60, a 0.25, B 5.7 T | 8 800 | 57 | 0.02 (floor) | 0 |
| chiral η 0.60 ε 0.70 n 4 (R 0.5, a 0.15) | 10 500 | 34 | 0.02 (floor) | 0.10–0.24 |
| chiral η 0.50 ε 0.55 n 3 (R 0.55, a 0.12) | 10 900 | 24 | 0.02 (floor) | 0.36–0.42 |

The criterion: the separation is shorted along the field if the connection length πR/ι
is shorter than the distance the shorting current — carried by the **electrons** — runs
while the column drifts one minor radius: ι ≳ (v_d,i/v_th,e)·πR/a, a few 10⁻³ at these
sizes, so any real transform (≥ 0.02, the floor used in the code) is enough. (The first
draft of this section used the ion thermal speed and quoted 0.11–0.26; the design agent
caught the inconsistency — see AGENT.md.) The criterion is satisfied by the chiral
configurations and not by the pair, which has ι = 0 exactly. (This rung is an estimate, not a simulation — a self-consistent electrostatic PIC
is the next thing to build here; the untwisted torus is its obvious validation case
because the answer is known.)

## 4. Rung 4 — forces and conductors

`ccsim/forces.py`, reactor preset (1 m ball, 364 kA per strand, conductor radius thinned
to 6.6 mm by the 1.5 cm clearance), untwisted 24 × 2 pair (`results/exp12_forces.json`,
fig15):

| quantity | value | against |
|---|---|---|
| B at the conductor surface | 17.5 T | magnetic pressure 121 MPa |
| peak force per length | 4.2 MN/m (p99 2.7, median 0.7) | μ₀I²/2πd at a crossing: 1.8 MN/m |
| contact/bending stress proxy f/2a | 320 MPa | annealed Cu 70, hard Cu 300, 316 SS 250, Inconel 1 000 MPa |
| current density | 2.7 × 10⁹ A/m² | steady Cu 2 × 10⁷ (×132), pulsed Cu 10⁹ (×2.6), NbTi 3 × 10⁸ (×9), Nb₃Sn 5 × 10⁸ (×5), REBCO 10⁹ (×2.7) |
| net force on each strand family | ±0.74 MN | the two families push apart; the frame must react it |

None of this is a showstopper; all of it says the reactor preset's 364 kA is 3–10× more
current than any conductor of that size carries. At 100 kA per strand (B_rms ≈ 1.6 T)
every ratio drops by 3.6 (forces by 13) and REBCO or Nb₃Sn work with margin; the physics
results above are scale-free in the sense that ι and the surfaces do not depend on the
current at all.

## 5. What "make sure it works" still needs

1. **A magnetic well.** Every chiral case is a hill. Add helical-axis excursion /
   triangularity on top of the modulation and re-run rung 1 — the optimiser loop is a
   few lines on `exp10`.
2. **Quasi-symmetry.** The 12 % collisionless loss in rung 2 will grow with collisions
   (1/ν regime). This is the standard stellarator optimisation problem; the mesh gives
   you two free functions on the torus (the two density modulations) to do it with.
3. **Finite β.** Vacuum surfaces are necessary, not sufficient; a 3D equilibrium
   (VMEC-class) is the next rung, and it will also give the real Mercier/ballooning
   answer to the hill.
4. **The self-consistent electrostatic test** of rung 3.
5. **Conductor engineering** at a current the materials allow, with the crossing loads
   reacted by a structure — the woven mesh's 1 150 crossings are 1 150 load points.

## 7. Session 8 — the corrected weave: surfaces versus the over/under offset δη

`experiments/exp13_weave_regimes.py` repeats the rung-1 Poincaré analysis (exact field, 24 seeds
from the magnetic axis, 40 toroidal turns) on the corrected geometry as a function of δη.

| η 0.60, ε 0.70, n 4, 24 circuits | weave | zones | clearance | surfaces | r_out | ι axis → edge | well |
|---|---|---|---|---|---|---|---|
| δη 0.005 | woven | 19 % | 0.0024 | 50 % | 0.168 | −0.10 → −0.25 | −0.23 |
| δη 0.010 | woven | 44 % | 0.0067 | 50 % | 0.166 | −0.10 → −0.25 | −0.22 |
| δη 0.015 | mixed | 68 % | 0.0100 | 46 % | 0.156 | −0.10 → −0.23 | −0.16 |
| δη 0.020 | mixed | 90 % | 0.0132 | 38 % | 0.132 | −0.10 → −0.18 | −0.13 |
| δη 0.030 | nested | 100 % | 0.0197 | 17 % | 0.072 | −0.19 → −0.22 | −0.13 |
| δη 0.050 | nested | 100 % | 0.0206 | 0 % | — | — | — |
| pure pair (ε 0), δη 0.010 | woven | 38 % | 0.0071 | 12 % | 0.022 | +0.01 → 0.00 | −0.06 |
| η 0.50, ε 0.55, n 3, δη 0.010 | woven | 43 % | 0.0082 | 50 % | 0.151 | -0.35 → -0.04 | -0.08 |
| η 0.50, ε 0.55, n 3, δη 0.030 | nested | 100 % | 0.0233 | 21 % | 0.069 | -0.42 → -0.23 | -0.06 |
| pure pair (ε 0), δη 0.030 | nested | 100 % | 0.0195 | 0 % | — | — | — |

(zones = the fraction of the strand where over/under is frozen because the other strand is too
close in the surface for a swap; r_out in ball units, tube radius 0.36.)

So: in the woven regime the chiral torus is exactly what §0 said (half the seeds on nested
surfaces out to 46 % of the tube radius, ι rising from 0.10 on axis to 0.25 at the edge, a
magnetic hill), and the pure pair has no transform (the achirality theorem is a statement about
the interleaved pair; a *nested* pair is not inversion-symmetric and does carry a coaxial poloidal
field in the gap between its shells — exp13 finds no surfaces there either). As δη grows the
interlock disappears turn by turn and the surfaces shrink with it; at δη 0.03 only a 0.07 core is
left and at 0.05 nothing. The over/under offset is therefore a design parameter with a hard
ceiling, and it is coupled to buildability: the true clearance of the woven torus is
2 δη |∂p/∂η| on the inboard side, 0.0067 at δη 0.01 (a 3 mm conductor at 1 m; the same-strand turn
spacing, 0.020, is the other ceiling). Rung 4 at 150 kA in a 3–4.5 mm conductor is
J ≈ 2.3–5 × 10⁹ A/m², two to five times the REBCO engineering limit — the meshed torus as a
single-conductor object is buildable at 2 T only with thinner tubes or fewer, fatter turns, and
that is now the agents' problem to solve on an honest geometry (AGENT.md §7–9: the woven δη 0.015
torus scores 17 with buildability 0.42; the same design at 1 T would score its full τ_c ≈ 41).

## 6. Files

`ccsim/fastfield.py`, `surfaces.py`, `guiding_center.py`, `polarization.py`, `forces.py`;
`geometry.hopf_mirror_pair` (tau0, mirror_rotation, deform, rank weave, closed knots),
`hopf_helical_pair`, `_modulated_strand`, `torus_deformation`, `twisted_mirror_pair`,
`hopf_torus_frame`; `experiments/exp10_surfaces.py` (parts 0/1, chiral, helical, best),
`exp11_transport.py`, `exp12_rung_figures.py`; figures fig12–fig15; the Lab's "Meshed torus II"
family, Poincaré button and Ledger tables. numba is now an optional dependency (the numpy
kernel is used when it is missing).
