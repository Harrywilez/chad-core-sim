# Leaderboard (score = τ_c × buildability, grid 25³, 160 markers, B_rms = 2 T)

`hbld` / `hscore` = *self-consistent* buildability and score: the same conductor radius
`a = min(1 cm, 0.45·clearance)` used in **both** buildability terms instead of the declared
`wire_radius` in one and `a` in the other (see notes.md F1). Ranking by `hscore` does not
change the top five.

| # | score | τ_c | build | hbld | hscore | clear. | J/J_REBCO | S4 | S20 | S120 | closed | design file |
|---|-------|-----|-------|------|--------|--------|-----------|----|-----|------|--------|-------------|
| 1 | **309.6** | 309.7 | 1.00 | 0.74 | 229.4 | 0.0211 | 0.46 | 0.38 | 0.24 | 0.21 | 0.00 | `designs/agent_knot_e0.88_N32.json` |
| 2 | 300.2 | 300.5 | 1.00 | 0.74 | 222.6 | 0.0177 | 0.52 | 0.39 | 0.23 | 0.21 | 0.00 | `designs/agent_knot_e0.88_N40.json` |
| 3 | 275.0 | 275.1 | 1.00 | 0.74 | 203.8 | 0.0211 | 0.40 | 0.34 | 0.23 | 0.19 | 0.00 | `designs/agent_knot_e0.85_N36.json` |
| 4 | 219.6 | 296.8 | 0.74 | 0.55 | 162.8 | 0.0155 | 1.35 | 0.69 | 0.31 | 0.20 | 0.00 | `designs/agent_knot_e1.15_N20.json` |
| 5 | 188.9 | 189.0 | 1.00 | 0.74 | 140.0 | 0.0219 | 0.86 | 0.62 | 0.36 | 0.21 | 0.00 | `designs/agent_knot_e1.12_N16.json` |
| 6 | 186.8 | 186.8 | 1.00 | 0.74 | 138.4 | 0.0170 | 0.94 | 0.56 | 0.33 | 0.21 | 0.00 | `designs/agent_knot_e1.05_N24.json` |
| 7 | 185.7 | 185.8 | 1.00 | 0.76 | 141.7 | 0.0229 | 0.85 | 0.57 | 0.33 | 0.21 | 0.00 | `designs/agent_knot_e1.10_N16.json` |
| 8 | 178.4 | 178.4 | 1.00 | 0.74 | 132.1 | 0.0195 | 0.87 | 0.58 | 0.36 | 0.21 | 0.00 | `designs/agent_knot_e1.10_N20.json` |
| 9 | 173.0 | 198.9 | 0.87 | 0.64 | 128.2 | 0.0190 | 1.15 | 0.67 | **0.43** | **0.26** | 0.00 | `designs/agent_knot_e1.15_N16.json` |
| 10 | 165.9 | 165.9 | 1.00 | 0.85 | 141.1 | 0.0255 | 0.41 | 0.35 | 0.18 | 0.15 | 0.00 | `designs/agent_knot_e0.80_N32.json` |
| 11 | 159.6 | 159.6 | 1.00 | 0.74 | 118.2 | 0.0144 | 0.99 | 0.63 | 0.34 | 0.26 | 0.00 | `designs/agent_knot_e1.00_N32.json` |
| 12 | 155.2 | 155.2 | 1.00 | 0.81 | 126.2 | 0.0244 | 0.47 | 0.36 | 0.21 | 0.14 | 0.00 | `designs/agent_knot_e0.85_N28.json` |
| 13 | 150.2 | 150.2 | 1.00 | 0.74 | 111.3 | 0.0185 | 0.59 | 0.46 | 0.31 | 0.23 | 0.00 | `designs/agent_knot_e0.92_N32.json` |
| 14 | 139.5 | 139.6 | 1.00 | 0.77 | 107.6 | 0.0231 | 0.41 | 0.40 | 0.24 | 0.18 | 0.00 | `designs/agent_knot_e0.85_N32.json` |
| 15 | 130.9 | 130.9 | 1.00 | 0.74 | 96.9 | 0.0163 | 0.77 | 0.48 | 0.29 | 0.21 | 0.00 | `designs/agent_knot_e0.96_N32.json` |
| 16 | 109.5 | 144.9 | 0.76 | 0.56 | 81.1 | 0.0124 | 1.32 | 0.60 | 0.39 | **0.27** | 0.00 | `designs/agent_knot_e1.05_N32.json` |
| 17 | 104.3 | 148.7 | 0.70 | 0.52 | 77.2 | 0.0139 | 1.43 | 0.67 | 0.36 | 0.25 | 0.00 | `designs/agent_knot_e1.10_N24.json` |
| 18 | 102.7 | 102.7 | 1.00 | 0.87 | 89.4 | 0.0261 | 0.41 | 0.28 | 0.15 | 0.11 | 0.00 | `designs/agent_knot_e0.75_N32.json` |
| 19 | **62.7** | 67.3 | 0.93 | 0.69 | 46.4 | 0.0148 | 1.07 | 0.46 | 0.26 | 0.11 | **0.51** | `designs/agent_m24w.json` |
| 20 | 60.5 | 60.5 | 1.00 | 0.74 | 44.8 | 0.0221 | 0.32 | 0.27 | 0.14 | 0.08 | 0.00 | `designs/agent_knot_eta0p70_N40_rev2.json` |
| 21 | 58.7 | 58.7 | 1.00 | 0.91 | 53.4 | 0.0273 | 0.39 | 0.28 | 0.11 | 0.08 | 0.00 | `designs/agent_knot32.json` |

Reference points from the pre-existing ledger (same score definition):

| score | τ_c | build | design |
|-------|-----|-------|--------|
| 24.1 | 91.4 | 0.264 | meshed torus II (hopfmesh η0.60 ε0.70 n4, wire 0.004) — previous best |
| 17.7 | 17.7 | 1.00 | hopf mirror pair (untwisted) |
| 5.3 | 5.3 | 1.00 | reference solenoid |
| 0.68 | 138.2 | 0.005 | precess × pair-axial-same (buildability ≈ 0) |

**Row 19 (`agent_m24w`) is the design I would actually recommend building** — see REPORT.md:
it is the only entry above 60 whose field has closed field lines (51 %), a 60-length
connection length and (from RUNGS §1, same family) verified nested flux surfaces with
ι ≈ 0.1–0.24. Every knot entry has `closed = 0.00`, connection length 2–5, mirror ratio
≈ 1.1–1.35 and, on the one `--surfaces` run made, **0 % flux surfaces**.
