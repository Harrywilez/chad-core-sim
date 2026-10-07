"""ccsim — Chad Core Simulator.

One SI-consistent stack that joins the three strands of the Chad Core Prime
work:

* geometry   — closed-wire windings: the recursive ∫dl → ∫∫dl → ∫∫∫dl concept
               from the notebook sketches, the three Codex realizations
               (phase-slip, precessing, Hopf drift), and reference coils.
* fields     — exact magnetostatic Maxwell solution for those windings
               (finite-segment Biot–Savart with finite conductor radius),
               the vector potential A, Faraday's induced E = −∂A/∂t for
               current ramps, and numerical ∇·B / ∇×B / Ampère checks.
* particles  — full-orbit Boris integration in SI with adaptive sub-stepping,
               wall/wire losses, and adiabatic-invariant diagnostics.
* collisions — test-particle Coulomb scattering against a Maxwellian background.
* fusion     — Bosch–Hale cross-sections and reactivities, Monte Carlo
               beam–target fusion events for tracked ions.
* chemistry  — reaction-network engine (mass action + Gillespie, Arrhenius with
               Eckart tunnelling), electron-impact rate coefficients from an
               EEDF, reduced CO₂/H₂ plasma network, Q-Surface import, and the
               coupling of ion wall flux to a surface network.
* presets    — benchtop, reactor-scale, and dimensionless configurations.

Everything in the physics modules states its own validity envelope; see
README.md for the fidelity ladder and what is *not* modelled.
"""

__version__ = "0.1.0"
