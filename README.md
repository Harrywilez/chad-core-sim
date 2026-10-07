# Chad Core Sim (`ccsim`)

One SI-consistent simulation stack that puts the three strands of the Chad Core
Prime work into a single folder and makes them talk to each other:

| strand | where it came from | what `ccsim` does with it |
|---|---|---|
| the coil concept (∫dl → ∫∫dl → ∫∫∫dl, "wire at 1×, 10×, 100×") | the three notebook pages | `geometry.recursive_winding` builds it to any depth with a rotation-minimising frame |
| the Codex closed-wire trials (phase-slip / precessing / Hopf drift, 6 Hopf-drift variants) | `~/.codex/visualizations/2026/08/31/…` | ported verbatim as generators; normalisation reproduced so results are directly comparable |
| the Fusion "Chat Core" lab (doublet equilibrium, Boris probes, Bosch–Hale D–T ledger) | `~/Documents/ChatGPT/Fusion` | Bosch–Hale reactivities cross-checked against it (2.7399e-22 m³/s at 15 keV matches to 5 digits); the ledger is now driven by *simulated* confinement instead of an assumed equilibrium |
| Q-Surface (reduced-order reaction chemistry) | `~/Documents/ChatGPT/Chemical reaction simulation` | its networks import directly; every `sech2` barrier becomes an exact Eckart (V0, ΔE, ω‡) barrier; a wall-flux map from particle tracing feeds the surface network |

## What "electromagnetically accurate" means here

`fields.py` is the exact magnetostatic solution of Maxwell's equations for the
windings: the Hanson–Hirshman finite-segment Biot–Savart formula (no softening
parameter outside the conductor, uniform-current interior inside it), the exact
vector potential **A** of the same filaments so that **B = ∇×A** holds, and the
induced electric field of a current ramp **E = −∂A/∂t**, which satisfies Faraday's
law identically. Every one of those statements is checked numerically in
`experiments/exp01_verify.py`:

* loop / solenoid / straight-wire / conductor-interior closed forms: ≤ 5e-5 relative
* ∇·B, ∇×B (outside conductors), ∇×A − B on real Chad-core windings: ≤ 7e-4 (finite-difference limited)
* Ampère loop integral around a conductor: exact to 2e-16
* Faraday EMF around a test loop vs −dΦ/dt: 3e-4

What is **not** included, deliberately: displacement current (irrelevant below
GHz at these sizes), the plasma's own currents and space charge (this is a
vacuum-field + test-particle stack), and any self-consistent MHD. An
electrostatic bias can be superposed from point charges / a uniform field.
The fidelity ladder therefore reads: exact vacuum EM → full-orbit test
particles with Coulomb collisions → Monte Carlo fusion events → reaction
networks. Each rung states its own validity envelope in its module docstring.

## Modules

```
ccsim/
  constants.py   CODATA constants, Species table
  geometry.py    windings: recursive (∫dl…), Codex ports, Hopf-drift variants, reference coils, remote return
  fields.py      exact Biot–Savart, A, induced E, electrostatics, FieldGrid cache, Maxwell checks, analytic references
  fieldlines.py  field-line tracing, mirror ratios, adiabatic trapped-fraction prediction
  particles.py   relativistic Boris pusher with adaptive sub-stepping, losses, invariants, ScaledGridField
  collisions.py  NRL test-particle Coulomb operator (Langevin MC), neutral-collision placeholder
  fusion.py      Bosch–Hale σ(E), ⟨σv⟩(T), beam–target Monte Carlo fusion events
  eckart.py      exact Eckart transmission, Γ(T), crossover temperature, (w, m) → ω‡
  chemistry.py   reaction-network engine (mass action + Gillespie), EEDF electron-impact rates,
                 CO₂/H₂ plasma network, illustrative Sabatier surface network, Q-Surface import, wall-flux map
  presets.py     benchtop / reactor / dimensionless scale presets
  runner.py      grids for the catalogue, per-preset particle screens, fusion ledger hooks
  surfaces.py    rung 1: batched field-line tracing on the exact field, Poincaré sections, magnetic axis, ι, well
  guiding_center.py  rung 2: guiding-centre orbits with pitch-angle scattering
  polarization.py    rung 3: polarization loss and the shorting criterion
  forces.py      rung 4: strand forces, stresses, conductor limits
  fastfield.py / fastorbits.py   numba kernels (Biot–Savart, Boris) — 20–25× faster, numpy fallbacks
  design.py      declarative designs: 16 component families (incl. the `cone` "tornado" coil) + transforms → MultiWinding
  evaluate.py    the grading function (score = τ_c × buildability) and the ledger (AGENT.md)
  export3d.py    3-D file export: OBJ polylines, CSV, watertight swept-tube meshes (binary STL, OBJ+MTL)
experiments/     exp01 verify · exp02 screen · exp03 topology · exp04 theory tests · exp05 chemistry · exp06 figures
                 exp07 continuation/assemblies · exp08 playbook · exp09 follow-ups · exp10–12 rungs · exp13 weave regimes
                 export_hopf_family*.py → exports/hopf_family/ (the drift winding and both mirror-pair versions as
                 equations + 3-D files + a note; README.md there)
agent/           PROTOCOL.md (for design agents) · optimize.py · rescore.py · notes, reports, leaderboard
designs/         design JSONs (seeds, agent designs, tornado plugs)
ui/              the Chad Core Lab (build.py, engine/app/body/head parts, serve.py live server, feed_export.py)
reference/       Codex factorial data and the exported Q-Surface network used for cross-checks
cache/           unit-scale field grids (B, A, wire distance), one per winding
results/         JSON outputs of every experiment · design_ledger.jsonl (every evaluation) · figures/  PNGs
```

## Run

```bash
python3 experiments/exp01_verify.py          # ~1 min: every module vs an independent reference
python3 experiments/build_grids.py           # ~6 min once: 13 unit-scale grids (cached)
python3 experiments/exp02_screen.py all      # ~1 h: particles in all configurations at all three scales
python3 experiments/exp03_topology.py        # ~6 min: field-line topology + adiabatic prediction
python3 experiments/exp04_theory_tests.py    # ~10 min: prediction vs orbits, adiabaticity sweep, Codex comparison, fusion ledger
python3 experiments/exp05_chemistry.py       # ~1 min: wall-flux maps → gas/surface networks, Q-Surface import
python3 experiments/exp06_figures.py         # figures for REPORT.md
python3 experiments/exp07_continuation_and_assemblies.py all   # hours: Hopf continuation + multi-core assemblies (33³)
python3 experiments/exp08_playbook.py 0 & python3 experiments/exp08_playbook.py 1 &   # ~2.5 h on two cores: the parameter sweep (PLAYBOOK.md)
python3 experiments/exp08_collect.py --print # merge the sweep into results/exp08_playbook.json and the Lab's ledger
python3 experiments/exp09_followups.py all   # ~1.5 h: energy scan, coarse-mesh mirror pairs, 60-transit run
python3 experiments/exp09_figures.py         # fig09–fig11
python3 ui/build.py                          # rebuild ui/chad_core_lab.html (+ standalone) from the parts
python3 experiments/exp10_surfaces.py best0 & python3 experiments/exp10_surfaces.py best1 &   # rung 1: flux surfaces / transform (exact field, numba)
python3 experiments/exp11_transport.py 0     # rung 2: guiding-centre confinement runs (index 0–5, ~30 min each)
python3 experiments/exp12_rung_figures.py    # fig12–fig15 + results/exp12_forces.json (rung 4)
python3 -m ccsim.evaluate designs/<design>.json [--surfaces]   # the grading function (AGENT.md); ~10–30 s
python3 agent/optimize.py designs/<design>.json --iters 20      # baseline optimiser; agent protocol in agent/PROTOCOL.md
python3 agent/rescore.py --surfaces 3        # re-score every design with the current geometry/score → agent/leaderboard.md (~40 min)
python3 experiments/exp13_weave_regimes.py   # rung 1 with the corrected (zone) weave: surfaces / ι vs the over/under offset δη
python3 ui/serve.py                          # the Lab with a LIVE Designs tab: polls results/design_ledger.jsonl while an agent runs
```

**Watching the agent live.** `python3 ui/serve.py` opens the Lab at http://localhost:8765/ui/chad_core_lab_standalone.html
and the Designs tab polls the ledger every 3 s: every `ccsim.evaluate` call (yours, `agent/optimize.py`'s, or a
Claude/Codex agent's) appears as a new row, the score chart extends, and with "follow the newest evaluation"
ticked the 3-D view rebuilds to each new design as it lands. Click any row to rebuild that design (same generators
and transforms as `ccsim.design`, cross-checked to 1e-13), then compute its field, lines, particles or Poincaré
section in the browser. The published Lab (claude.ai artifact) has the same tab fed from a shared live feed that
a Claude session pushes evaluations into (`ui/feed_export.py` prepares the rows).

See **AGENT.md** for the design-loop (components → assemblies → score → agent), what two
agent runs found, and the session-8 geometry correction (the weave is now a *zone weave* and the
clearance is exact — the earlier "open weave" scores were an artifact), **PLAYBOOK.md** for what the sweep found (the woven Hopf mirror pair — the "meshed torus" — is
the standout of the single-particle screen), **RUNGS.md** for whether it would actually work
(it needs a chiral weave to have flux surfaces; that version is the Lab's "Meshed torus II"
family), and **ui/chad_core_lab_standalone.html** for the browser Lab.
Optional: `pip install numba` makes the exact-field tracing 25× faster (rungs 1–2).

Dependencies: NumPy and SciPy only (matplotlib for figures). Python ≥ 3.9.

## Using it as a library

```python
from ccsim.geometry import recursive_winding, close_with_return, Winding
from ccsim.fields import FieldGrid, maxwell_checks
from ccsim.particles import ScaledGridField, monoenergetic_ensemble, run_orbits
from ccsim.constants import DEUTERON

pts = recursive_winding("circle", base_radius_m=0.72, levels=[(0.24, 12), (0.30, 6)])
closed, n_active = close_with_return(pts, 3.8)
w = Winding("my coil", closed, current_A=1.0, wire_radius_m=0.005, closed=True, active_count=n_active)
grid = FieldGrid(w, half_extent_m=0.87, n=41, with_A=True)        # unit scale, once
field = ScaledGridField(grid, scale=0.18, current_A=1000.0)        # benchtop: 18 cm, 1 kA
ens = monoenergetic_ensemble(DEUTERON, 256, 10.0, positions, rng)
res = run_orbits(ens, field, t_max=2e-4, dt=1e-7, wall_radius=0.148)
```

## Fidelity statements you should read before quoting a number

* **Retention is a finite-window statistic**, exactly as in the Codex trials. It becomes physically meaningful only when compared with the adiabatic prediction from field-line topology (`fieldlines.py`) — see REPORT.md, tests T1/T2.
* **Fusion rates** are beam–target Monte Carlo on a fixed Maxwellian background: correct physics for tracked ions, but no self-heating, no alpha slowing-down feedback, no transport solution — so no plasma Q.
* **Chemistry**: the O + H₂ and OH + H₂ rates are literature values; three-body and electron-impact dissociation cross-sections are labelled PLACEHOLDER shapes; the Sabatier surface barriers are ILLUSTRATIVE and are the numbers your MACE CI-NEB chain is meant to supply as (V0, ΔE, ω‡).
* **Conductor thinning**: when a winding's clearance is smaller than the preset conductor diameter, the conductor is thinned to 45 % of the clearance and the resulting Joule power / current density are reported — several Codex geometries become physically unbuildable at benchtop scale for this reason, and the report says so.

## Companion repository

The multi-agent design search that used this simulator — the Together room where Codex, Claude and Astra
agents combined components and hunted for confinement, the Crucible rig-study room, the direct-vs-construction
pilot, the fable design search (413 evaluations) and its audit — is archived in the companion repository
`together` (`archive/README.md` there lists the rooms; every room's events, messages and files are exported as
plain text). The 413 fable-search evaluations are also appended to this repository's `results/design_ledger.jsonl`.
