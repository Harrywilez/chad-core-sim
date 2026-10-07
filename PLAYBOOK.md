# Chad Core Lab — playbook (session 5, 2026-09-02)

> **Read RUNGS.md after this.** The single-particle screen below made the meshed torus the
> winner; the flux-surface rung (session 6) shows the pure mirror pair has no rotational
> transform and would not hold a plasma — and how a chiral weave fixes that.

*What happened when I played with the Lab: a parameter sweep over every family, three new
geometries, and the Hopf-torus "meshed torus" mode you asked for. Everything here is
reproducible from `experiments/exp08_playbook.py`, `exp09_followups.py`, `exp09_figures.py`;
the numbers live in `results/exp08_playbook.json`, `results/exp09_*.json` and in the Lab's
Ledger tab ("Playbook sweep").*

## 0. The one-paragraph answer

The best thing in the box is the new one: **the Hopf torus plus its mirror image, woven
together, with the two strands carrying opposite currents.** At reactor scale (1 m ball,
364 kA per strand) it keeps **46 % of 15 keV deuterons after 20 wall transits** (51 % with a
fatter torus, η = 0.80), against 24–29 % for the best precessing-loop pair and ≤ 19 % for any
single Hopf winding. More importantly the survival curve goes *flat*: 5 keV → 47 %, 15 keV →
46 %, 45 keV → 42 %. The mirror-type configurations lose retention as the Larmor radius grows
(the precessing pair drops 42 → 26 → 10 % over the same energies); the single Hopf windings are
also flat in energy but flat at a low level (hopf continued 720°/48: 18 → 14 → 12 %), because
their loss is plain streaming along long open lines, which in units of wall transits does not
depend on speed. The meshed torus is different in kind: the mirror strand cancels the toroidal
current, leaving a poloidal current sheet — a **toroidal solenoid** — so the field lines inside
the tube close on themselves and there is no loss cone. Whatever starts inside the tube stays
(15–25 % of field-line seeds are on closed lines; the rest of the survivors sit on lines that
are long enough to matter). What starts outside is lost in about one transit. It still leaks
— a toroidal field has no rotational transform, so the vertical drift takes the tube population
out on a ~60-transit e-fold (S60 = 26 %) — but that is 5–7× slower than any mirror here.

## 1. The meshed-torus mode (what you asked for)

*Family → "Hopf torus + mirror image (meshed torus)" in the Lab.*

The mirror image (z → −z) of a Hopf fibre is an **anti-Hopf fibre**: it lies on the same torus
but winds the other way (a (1,−1) curve instead of (1,1)). So the two strands intersect each
other — with N circuits each they cross 2·(N² − r²/4) times (1152 crossings for 24 × 1) and
make a rhombic mesh on the torus. Two ways to fit them together:

* **woven** — at every crossing one strand is pushed outward (η + δη) and the other inward
  (η − δη), alternating along both strands like a plain weave. The crossings are found
  analytically: with p = N + r/2, q = N − r/2 they sit at u = (j/p + k/q)/2 on the original and
  v = (k/q − j/p)/2 on the mirror, and the over/under pattern is (−1)^k (98 % alternation
  along both strands). Mutual clearance 0.0116 ball units at δη = 0.025 (1.2 cm at 1 m; the
  conductor radius is 0.5 cm).
* **nested** — the original on the torus η + δη, the mirror on η − δη (clearance 0.0137).

*(Reproducibility note, added in session 6: this sweep used strand phase τ₀ = π/2, open
strands with remote returns, and the (−1)^k weave with independent profiles on the two
strands — `hopf_mirror_pair(..., tau0=np.pi/2, mirror_rotation=0.0)` in the current code; the
Lab's "Hopf torus + mirror image" family is exactly that. The defaults are now τ₀ = 0 and
mirror rotation π, and RUNGS.md explains why.)*

And two ways to feed them:

* **opposed** (mirror carries −I): toroidal currents cancel, poloidal currents add → toroidal
  field inside the tube, almost nothing outside (fig11, middle). This is the confining one.
* **same**: poloidal currents cancel → a fat current loop, dipole-like field through the hole
  (fig11, right). A mirror, S20 ≈ 19 %.

| meshed torus (η, circuits, fit, feed) | closed lines | S4 | S8 | S20 | B_rms | conductor clearance |
|---|---|---|---|---|---|---|
| 0.80, 24, woven, opposed | 23 % | 64 | 57 | **51 %** | 6.5 T | 0.010 |
| 0.70, 24, woven, opposed | 15 % | 54 | 48 | **46 %** | 5.7 T | 0.014 |
| 0.70, 36, woven, opposed | 25 % | 53 | 49 | 44 % | 8.6 T | 0.012 |
| 0.70, 24, nested, opposed | 2 % | 52 | 47 | 41 % | 5.7 T | 0.014 |
| 0.70, 16, woven, opposed | 11 % | 49 | 44 | 36 % | 3.8 T | 0.014 |
| 0.60, 24, woven, opposed | 15 % | 43 | 37 | 32 % | 4.8 T | 0.016 |
| 0.80, 24, woven, same | 0 % | 42 | 36 | 24 % | 5.4 T | 0.010 |
| 0.70, 24, woven, same | 0 % | 31 | 28 | 19 % | 6.8 T | 0.014 |
| 0.70, 24, mirror at −0.5 I | 0 % | 44 | 37 | 21 % | | |
| 0.70, 24, mirror at 0 (single strand) | 0 % | 36 | 23 | 11 % | | |
| 0.70, 8 / 6 / 4 circuits, opposed | 4 / 3 / 0 % | 46 / 44 / 41 | | 13 / 8 / 4 % | | |

What the table says:

* **Fatter torus is better** (η 0.60 → 0.70 → 0.80: 32 → 46 → 51 %) because the tube takes up
  more of the ball, so more of the starting population is inside it. η = 0.80 is at the
  clearance limit for a 0.5 cm conductor; η ≈ 0.85 would need a thinner wire.
* **Finer mesh is better** up to a point (16 → 24 → 36 circuits: 36 → 46 → 44 %). The coarse
  meshes (4–8 circuits) are not solenoids at all — the field ripples out between strands and
  the lines open (closed fraction ≤ 4 %). Rotational transform measured on the long lines of
  the 4-circuit pair is ι ≈ 0.2 per toroidal turn, but with no closed surfaces it doesn't help.
* **Woven beats nested** (46 vs 41 %) — small but consistent with the woven pair being a
  cleaner current sheet (the two strands sit on the *same* mean torus instead of two).
* **Unequal currents don't help.** Feeding the mirror at −0.5 I leaves a net toroidal current
  whose poloidal field opens the lines (closed 0 %, S20 21 %). This is the vacuum version of
  the tokamak problem: a smooth double current sheet gives you a toroidal field *or* a
  poloidal field inside the tube, never a rotational transform — that needs a discrete helical
  ripple (stellarator) or plasma current.
* **The energy scan is the signature.** S20 = 47 / 46 / 42 % at 5 / 15 / 45 keV (ρ_L = 2.5 /
  4.4 / 7.6 mm at 1 m). Only the 45 keV point shows the ∇B/curvature drift beginning to erode
  the tube population. Median loss time of the lost particles is ~1.1 transits at every energy:
  the loss is a prompt topological split, not a slow leak. (Reading the scan in general: a
  loss that is flat in *transits* across energies is streaming or topology; a loss that grows
  with energy is drift. The precess pair is drift, the single Hopf strand is streaming, the
  meshed torus is topology with a little drift on top.)

* **It is not a perfect trap.** The 60-transit run (η 0.70, 24, woven, opposed, 15 keV) gives
  S20 / S40 / S60 = 46 / 36 / 26 %. A pure toroidal field has no rotational transform, so the
  vertical ∇B/curvature drift still walks the tube population out — but with an e-folding of
  roughly 60–70 transits instead of the ~10 of the best mirror (precess pair: 69 → 26 % between
  4 and 20 transits). Fixing that is the whole history of toroidal confinement: you need
  either a plasma current (tokamak) or a helical ripple with closed flux surfaces
  (stellarator). The coarse-mesh experiments in the table are the first crude step toward
  the second option and they did not find surfaces; that is the obvious next thing to try
  properly (few strands with an l = 2 helical pitch on top of the mesh).

## 2. Precessing-loop pairs (the mirror family's best)

| precess pair variant | Rmir p90 | pred | S4 | S8 | S20 | core conc. |
|---|---|---|---|---|---|---|
| pair-axial-same, core 0.55 | 12.9 | 0.68 | 59 | 50 | **36 %** | 0.39 |
| pair-side-crossed, core 0.55 | 7.8 | 0.56 | 62 | 51 | 29 % | 0.28 |
| pair-side-crossed, rot 12° | 15.0 | 0.66 | 64 | 53 | 28 % | 0.31 |
| pair-side-crossed, core 0.45 (baseline) | 16.7 | 0.65 | 69 | 57 | 26 % | 0.29 |
| pair-side-crossed, inward 0.15 | 10.1 | 0.59 | 66 | 57 | 26 % | 0.28 |
| pair-side-crossed, inward 0.30 | 25.3 | 0.65 | 68 | 54 | 21 % | 0.32 |
| pair-side-crossed, rot 36° | 13.1 | 0.64 | **72** | 53 | 21 % | 0.25 |
| pair-side-crossed, core 0.35 | 32.2 | 0.71 | 65 | 43 | 12 % | 0.29 |

* These are the best *mirrors*: S4 up to 72 %, in line with the adiabatic prediction
  (√(1 − B_seed/B_mirror)), and then the drift tax: S20 is 2–5× lower than S4.
* **Bigger cores hold longer** (core scale 0.35 → 0.45 → 0.55: S20 12 → 26 → 29 %) even though
  the mirror ratio drops (32 → 17 → 8), because the drift-loss time scales with L_B/ρ_L and a
  bigger core has a longer field-scale length. Same story in the energy scan: 42 / 26 / 10 %
  at 5 / 15 / 45 keV — textbook ρ_L scaling (fig10).
* **Rotation 12° vs 36°**: slower precession gives the higher plateau (28 %), faster gives
  the higher initial capture (72 % at S4) — the fast-precessing core has more mirror throats
  but they are leakier.
* The precessing winding's own clearance is its Achilles heel: 0.0011–0.0023 ball units, i.e.
  a conductor **1–2 mm thick at 1 m** carrying 364 kA. The Joule ledger says 0.07–1.3 TW of
  copper losses at reactor scale, against ~4 GW for the meshed torus (1.4 cm wire). Neither is
  a real magnet design, but the ratio is the point: the precessing pair is a great field for
  a hard geometry.
* **pair-axial-same at core 0.55** (36 %) has the highest core concentration of the mirrors
  (39 % of the live particles inside r < 0.35 on average) — if "convergence" means pulling the
  population into the centre, this is the pick among the mirrors.

## 3. The single Hopf family (drift / continued / torus)

| configuration | Rmir p90 | pred | S4 | S8 | S20 |
|---|---|---|---|---|---|
| hopf torus η = 0.78, 24 × 1 | 1.38 | 0.23 | 38 | 29 | 19 % |
| hopf torus η = 0.70, 25 × 3 | 1.50 | 0.24 | 31 | 24 | 15 % |
| hopf continued 360° / 24 | 1.29 | 0.19 | 38 | 28 | 15 % |
| hopf continued 720° / 48 (η→0.72) | 1.16 | 0.18 | 28 | 21 | 14 % |
| hopf continued 720° / 48 | 1.19 | 0.16 | 25 | 20 | 14 % |
| hopf torus η = 0.70, 24 × 1 | 1.37 | 0.22 | 32 | 21 | 12 % |
| hopf continued 720° / 32 | 1.15 | 0.16 | 22 | 18 | 11 % |
| hopf torus η = 0.62 / 0.55 | 1.3 | 0.21 | 24 / 20 | | 11 / 7 % |
| hopf continued 540° / 36 | 1.54 | 0.14 | 22 | 13 | 9 % |
| hopf continued 720° / 64 | 1.12 | 0.13 | 20 | 12 | 9 % |
| hopf continued 720° / 48 (η→0.56) | 1.12 | 0.12 | 16 | 12 | 4 % |
| hopf drift Shallow 180°, 16 / 10 / 6 circuits | 4.6–5.4 | 0.23 | 22 / 21 / 17 | | 5 / 4 / 0 % |

* A single Hopf strand is a (N,N) torus curve — 45° pitch — so it makes a toroidal *and* a
  poloidal field of similar size, and its mirror ratios are barely above 1. It never
  confines well on its own (S20 ≤ 19 %); it is the raw material for the meshed torus.
* **More circuits is not better** past ~24–48: 720° / 32 → 48 → 64 gives 11 → 14 → 9 %. The
  field gets stronger (5.6 → 8.4 → 11.2 T) but the geometry gets *smoother*, so the mirror
  ratio falls (1.15 → 1.19 → 1.12). Retention here is set by mirror ratio, not field strength
  (at fixed ρ_L/L_B the Boris orbits don't care about |B|).
* Letting the continuation end deeper (η → 0.56) is bad (4 %); ending shallower (η → 0.72)
  slightly better than the default. The 360° / 24 sweep is as good as 720° / 48 with half the
  conductor.
* Fatter torus helps here too (η 0.55 → 0.78: 7 → 19 %), for the same reason as in §1 — more
  of the ball is inside the tube.
* The Codex-style drift winding (Shallow 180°) is the weakest of the family at every circuit
  count: high mirror ratio (4.6–5.4) but tiny field scale, so the drift empties it by 20
  transits.

## 4. New geometries

| configuration | Rmir p90 | pred | S4 | S8 | S20 | core conc. | B_rms |
|---|---|---|---|---|---|---|---|
| baseball × pair-side-crossed | 25.6 | **0.81** | 57 | 40 | 8 % | 0.36 | 1.7 T |
| sphere winding, 24 turns (uniform B) | 1.03 | 0.08 | 26 | 16 | 10 % | 0.08 | 5.1 T |
| yin-yang pair (4 + 4 turns) | 2.1 | 0.23 | 20 | 12 | 7 % | 0.03 | 2.3 T |
| picket fence, 7 rings | 4.7 | 0.62 | 19 | 8 | 2 % | **0.41** | 0.31 T |
| picket fence, 5 rings | 8.1 | 0.63 | 21 | 7 | 0 % | 0.22 | 0.37 T |
| baseball seam, 60° | 1.01 | 0.05 | 6 | 2 | 1 % | **0.68** | 2.4 T |
| baseball seam, 40° | 1.00 | 0.02 | 1 | 0 | 0 % | 0.09 | 2.3 T |
| Hopf link (two linked rings), same / opposed | 3.1 / 2.7 | 0.30 | 11 / 12 | 3 / 1 | 0 % | 0.08 | 0.54 T |

* **Baseball-seam coil**: on its own it is the worst thing in the catalogue — the seeds sit in
  the coil's own high-field region and every line runs downhill to the wall (mirror ratio
  1.00). Livermore's baseball coils worked because the plasma sat *inside* a coil much
  bigger than itself. As a **pair-side-crossed assembly** it has the highest adiabatic
  prediction of anything (0.81, mirror ratio 26) and S4 = 57 %, but the drift tax is
  brutal (S20 = 8 %). Good capture, poor hold.
* **Sphere winding**: the uniform-field check. Turns equally spaced in z give a surface
  current K ∝ sin θ, whose interior field is exactly uniform; with the winding cut off at
  |z| = 0.94 R the closed form is B = μ₀NI(c − c³/3)/(2cR), c = 0.94, and the Biot–Savart
  centre field matches it to 0.06 % (now part of `exp01_verify`; uniform to 2 % over the
  interior sample, B_max/B_rms = 1.0 in the ROI). Confinement
  is exactly what a uniform field gives: the 10 % that survive are the near-perpendicular
  pitch angles; nothing drifts because there is no gradient.
* **Picket fence**: the cusps cancel most of the field (B_rms 0.3 T for the same current) but
  it has the highest *core concentration* of any multi-circuit geometry (41 %) — a
  ring-cusp really does push what it holds into the middle. Survival is short.
* **Yin-yang** and **Hopf link**: pretty, not useful here. The link is a helicity toy — B·A
  is non-zero, but at this size the field is 0.5 T and mirror ratios ~3 with tiny L_B.

## 5. Rules of thumb that came out of this

1. **Closed lines beat mirror ratio.** Every mirror in the catalogue tops out at S4 ≈ 0.7 and
   loses 2–5× of that by 20 transits; the only configurations with a flat tail are the ones
   with closed field lines (the meshed torus). If the target is long confinement, cancel the
   toroidal current and keep the poloidal one.
2. **The drift tax scales with ρ_L/L_B.** Bigger core scale, lower energy, or a longer
   field-scale length all help mirrors by the same factor (precess core 0.35 → 0.55 doubled
   S20; 45 → 5 keV quadrupled it). Field *strength* on its own does nothing for a mirror.
3. **The Hopf family wants a mirror partner, not more circuits.** Past ~24 circuits per
   revolution the single strand just becomes a smoother, weaker mirror.
4. **Fatter tori win** (η ≈ 0.8 for the meshed torus; 0.78 for the single torus) because
   retention is really "fraction of the vessel that is inside the tube".
5. **Capture and hold are different metrics.** baseball-pair / precess rot 36° capture best
   (S4 ≈ 0.6–0.7); meshed torus holds best (S20 ≈ 0.5). Core concentration is a third axis:
   picket fence and the axial precess pair pull the population inward.

## 6. Caveats and the numbers I trust least

* 25³ grids (exp08/09) versus 33³ (exp07): the same configuration moves by a few points
  (hopf continued 720° / 48: S20 = 0.19 at 33³ / 256 particles vs 0.14 at 25³ / 160). Treat
  differences of < 5 points as noise; 160 particles is ±4 % (1σ) at S ≈ 0.5.
* No collisions, no plasma back-reaction, no wall recycling; single-energy D⁺. The Ledger's
  fusion column is not part of this sweep.
* The reactor preset feeds 364 kA into copper conductors sized by the geometry's clearance —
  the Joule numbers are physical for that copper but nobody would build it that way. Use them
  as a *relative* cost of the geometry's tightness.
* Rotational-transform measurement (exp09, `transform`) uses a round-torus estimate of the
  centre circle from the strand's ρ-range; it is fine for ι ≈ 0.2 and meaningless when the
  lines are open.
* The exp07 leftovers (precess triad-alternating / tetra / octa / cube; solenoid assemblies)
  were not finished — the background processes did not survive the container pause, and the
  meshed-torus result made them less interesting. `experiments/exp07_continuation_and_assemblies.py assemblies-precess`
  resumes them from the cache.

## 7. Where things are

* Lab (browser): family "Hopf torus + mirror image (meshed torus)" and the five new families
  (baseball, yin-yang, picket fence, sphere winding, Hopf link); Ledger tab → "Playbook sweep"
  lists every row above, click to load its exact parameters. The Ampère check now picks the
  most isolated conductor segment for its loop, so multi-circuit meshes pass it.
* Python: `ccsim.geometry.hopf_mirror_pair`, `hopf_mirror_windings`, `baseball_seam`,
  `yin_yang`, `picket_fence`, `sphere_winding`, `hopf_link`, `multi_windings`.
* Figures: `figures/fig09_playbook_ranking.png`, `fig10_energy_scan.png`,
  `fig11_meshed_torus_anatomy.png`.
