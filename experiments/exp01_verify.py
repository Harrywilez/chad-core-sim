"""Verification suite: every physics module against an independent reference.

  fields     segment Biot–Savart vs analytic loop / solenoid / wire / conductor interior;
             ∇·B, ∇×B, ∇×A−B, Ampère loop on a real Chad-core winding;
             Faraday: ∮E_ind·dl = −dΦ/dt for a current ramp.
  particles  gyroradius and gyrofrequency in a uniform field; μ conservation in a
             magnetic mirror (two-loop bottle) and the loss-cone prediction;
             ∇B drift speed vs the guiding-centre formula.
  collisions slowing-down of a 3.5 MeV alpha on a 15 keV electron background vs the
             analytic (small-x) NRL rate; energy relaxation of a test Maxwellian.
  fusion     Bosch–Hale σ vs ⟨σv⟩ fits, reference values, beam–target MC vs n σ v.
  chemistry  element/site balance, Arrhenius + Eckart Γ(T) against the Q-Surface numbers,
             EEDF integration against a closed form.
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from ccsim.chemistry import co2_h2_plasma_network, eedf_rate_coefficient, sabatier_surface_network
from ccsim.collisions import Background, CoulombOperator, relaxation_rates, slowing_down_time
from ccsim.constants import ALPHA, DEUTERON, ELECTRON, E_CHARGE, EPS0, M_E, MU0
from ccsim.eckart import crossover_temperature_K, imag_freq_from_width, tunnelling_factor
from ccsim.fields import FieldGrid, analytic_references, biot_savart, induced_electric_field, loop_field_analytic, maxwell_checks
from ccsim.fusion import FusionEventSampler, cross_section_m2, reactivity_m3_s, verify as fusion_verify
from ccsim.geometry import Winding, chadcore_catalogue, circular_loop
from ccsim.particles import Ensemble, ScaledGridField, run_orbits

out = {}
t_start = time.time()

# ---------------------------------------------------------------- fields
out["fields_analytic"] = analytic_references()
cat = chadcore_catalogue(1.0, 0.005, 1.0, include_hopf_drift=False)
w = cat[1]   # recursive L3
rng = np.random.default_rng(3)
pts = rng.uniform(-0.6, 0.6, size=(80, 3))
out["fields_maxwell_recursive_L3"] = maxwell_checks(w, pts)
w2 = cat[5]  # codex hopf
out["fields_maxwell_codex_hopf"] = maxwell_checks(w2, pts)

# Faraday: EMF around a test loop = -dΦ/dt.  Φ through a small loop of radius r0 at the origin, plane z=0.
loop = Winding("loop", circular_loop(0.3, n=1440), 1000.0, 1e-3)
dIdt = 1.0e6   # A/s
r0 = 0.05
th = np.linspace(0, 2 * np.pi, 721)
ring = np.column_stack((r0 * np.cos(th), r0 * np.sin(th), np.zeros_like(th)))
E = induced_electric_field(ring, loop, dIdt)
tang = np.column_stack((-np.sin(th), np.cos(th), np.zeros_like(th)))
emf = np.trapezoid(np.sum(E * tang, axis=1) * r0, th)
# flux through the ring: integrate Bz over the disc (polar quadrature) at unit current, times dI/dt
rr = np.linspace(0, r0, 60)[1:]
pp = np.linspace(0, 2 * np.pi, 241)
R, P = np.meshgrid(rr, pp, indexing="ij")
disc = np.column_stack((R.ravel() * np.cos(P.ravel()), R.ravel() * np.sin(P.ravel()), np.zeros(R.size)))
Bz = biot_savart(disc, loop, current_A=1.0)[:, 2].reshape(R.shape)
flux_per_amp = np.trapezoid(np.trapezoid(Bz * R, pp, axis=1), rr)
out["faraday_emf_rel_err"] = float(abs(emf / (-dIdt * flux_per_amp) - 1.0))

# ---------------------------------------------------------------- particles
class Uniform:
    def __init__(s, B):
        s.B0 = np.array(B); s.half_extent = 1e9; s.wire_radius = 0.0; s.scale = 1.0
        s.current_ramp = None; s.electrostatics = None
        class G: dx = 1e-3
        s.grid = G()
    def B_at(s, x, t=0): return np.tile(s.B0, (len(x), 1))
    def E_at(s, x, t=0): return np.zeros_like(x)
    def Bmag_at(s, x): return np.full(len(x), np.linalg.norm(s.B0))
    def wire_distance_at(s, x): return np.full(len(x), 1e9)
    def exact_wire_distance(s, x): return np.full(len(x), 1e9)

B0 = 1.0
f = Uniform([0, 0, B0])
vperp = 1e6
omega = E_CHARGE * B0 / DEUTERON.mass_kg
T_gyro = 2 * math.pi / omega
ens = Ensemble(DEUTERON, np.zeros((1, 3)), np.array([[vperp, 0.0, 0.0]]))
res = run_orbits(ens, f, t_max=20 * T_gyro, dt=T_gyro / 60, wall_radius=1.0, omega_dt_max=0.2, sample_every=1)
xy = np.array([s[1][0, :2] for s in res.extra["samples"]])
rho_num = 0.5 * (xy[:, 0].max() - xy[:, 0].min())
rho_th = DEUTERON.mass_kg * vperp / (E_CHARGE * B0)
# gyro-period from zero crossings of y
y = xy[:, 1] - xy[:, 1].mean()   # orbit is not centred on the origin
ts = np.array([s[0] for s in res.extra["samples"]])
crossings = np.flatnonzero(np.diff(np.sign(y)) != 0)
tc = ts[crossings] - y[crossings] * (ts[crossings + 1] - ts[crossings]) / (y[crossings + 1] - y[crossings])
period_num = 2 * (tc[-1] - tc[0]) / (len(tc) - 1)
out["particles_uniform"] = {"gyroradius_rel_err": float(abs(rho_num / rho_th - 1)), "gyroperiod_rel_err": float(abs(period_num / T_gyro - 1)),
                            "boris_phase_error_theory": float(abs(2 * math.atan(omega * T_gyro / 60 / 2) / (omega * T_gyro / 60) - 1)),
                            "energy_drift": res.energy_drift_rel_max, "mu_variation": res.mu_variation_median}

# magnetic mirror: two coaxial loops; predict trapping from the loss cone sin^2(theta_c) = B_min/B_max
from ccsim.geometry import circular_loop as _cl
Rm = 0.2
l1, l2 = _cl(Rm, 720, center=(0, 0, -0.3)), _cl(Rm, 720, center=(0, 0, 0.3))
bottle_pts = np.vstack((l1, l2, l1[:1]))   # up-and-back jumper cancels exactly
bottle = Winding("bottle", bottle_pts, 2.0e4, 2e-3)
grid = FieldGrid(bottle, 0.45, 61)
fld = ScaledGridField(grid, 1.0, 2.0e4)
Bmin = float(fld.Bmag_at(np.zeros((1, 3)))[0])
Bmax = float(fld.Bmag_at(np.array([[0, 0, 0.3]]))[0])
theta_c = math.asin(math.sqrt(Bmin / Bmax))
n = 400
rng = np.random.default_rng(5)
pitch = rng.uniform(0, math.pi, n)
speed = math.sqrt(2 * 10.0 * E_CHARGE / DEUTERON.mass_kg)   # 10 eV deuterons (rho_L ~ 2 cm << bottle)
phi = rng.uniform(0, 2 * math.pi, n)
v = speed * np.column_stack((np.sin(pitch) * np.cos(phi), np.sin(pitch) * np.sin(phi), np.cos(pitch)))
x0 = np.zeros((n, 3)) + rng.normal(0, 0.005, size=(n, 3))
ens = Ensemble(DEUTERON, x0, v)
res = run_orbits(ens, fld, t_max=30 * 0.6 / speed, dt=3e-8, wall_radius=0.44, omega_dt_max=0.2, wall_shape="cylinder", wall_half_length=0.42)
trapped_pred = (pitch > theta_c) & (pitch < math.pi - theta_c)
agree = np.mean(trapped_pred == ens.alive)
out["particles_mirror"] = {"B_min": Bmin, "B_max": Bmax, "loss_cone_deg": math.degrees(theta_c), "predicted_trapped_fraction": float(trapped_pred.mean()),
                           "observed_trapped_fraction": float(ens.alive.mean()), "per_particle_agreement": float(agree),
                           "mu_variation_median": res.mu_variation_median, "energy_drift": res.energy_drift_rel_max, "adiabaticity": res.adiabaticity_median}

# grad-B drift: field B = B0 (1 + x/L) z-hat  -> v_d = (m v_perp^2 / (2 q B)) * (1/L)
class GradB(Uniform):
    def __init__(s, B0, L): super().__init__([0, 0, B0]); s.L = L
    def B_at(s, x, t=0): return np.column_stack((np.zeros(len(x)), np.zeros(len(x)), s.B0[2] * (1 + x[:, 0] / s.L)))
    def Bmag_at(s, x): return s.B0[2] * (1 + x[:, 0] / s.L)
L = 5.0
g = GradB(1.0, L)
ens = Ensemble(DEUTERON, np.zeros((1, 3)), np.array([[vperp, 0.0, 0.0]]))
t_max = 400 * T_gyro
res = run_orbits(ens, g, t_max=t_max, dt=T_gyro / 80, wall_radius=10.0, omega_dt_max=0.1, sample_every=1)
ts = np.array([s[0] for s in res.extra["samples"]]); ys = np.array([s[1][0, 1] for s in res.extra["samples"]])
vd_num = np.polyfit(ts, ys, 1)[0]
vd_th = (DEUTERON.mass_kg * vperp**2 / (2 * E_CHARGE * 1.0)) / L   # v_d = -mu (grad B x B)/(q B^2) = +y for +q with B rising in +x
out["particles_gradB_drift"] = {"v_drift_numeric": float(vd_num), "v_drift_theory": float(vd_th), "rel_err": float(abs(vd_num / vd_th - 1))}

# ---------------------------------------------------------------- collisions
bg_e = Background(ELECTRON, 1e20, 15e3)
tau_s = slowing_down_time(ALPHA, 3.5e6, bg_e)
# small-x limit of the NRL rates: nu_s -> (1+m_a/m_e) * (4/(3 sqrt(pi))) x^{3/2} nu_0
v_a = math.sqrt(2 * 3.5e6 * E_CHARGE / ALPHA.mass_kg)
x = M_E * v_a**2 / (2 * 15e3 * E_CHARGE)
lnL = 17.5
nu0 = bg_e.density_m3 * (2 * E_CHARGE) ** 2 * E_CHARGE**2 * lnL / (4 * math.pi * EPS0**2 * ALPHA.mass_kg**2 * v_a**3)
nu_s_small = (1 + ALPHA.mass_kg / M_E) * (4 / (3 * math.sqrt(math.pi))) * x**1.5 * nu0
out["collisions"] = {"alpha_tau_s_s": tau_s, "small_x_limit_tau_s": 1 / nu_s_small, "rel_diff": abs(tau_s * nu_s_small - 1), "x": x}
# Monte Carlo energy relaxation: 1 keV deuterons on 1 keV deuteron background should keep <E> ~ (3/2)T... test mean energy drift over 0.2 tau
rng = np.random.default_rng(11)
bg_i = Background(DEUTERON, 1e20, 1e3)
n = 2000
ens = Ensemble(DEUTERON, np.zeros((n, 3)), rng.normal(0, math.sqrt(1e3 * E_CHARGE / DEUTERON.mass_kg), size=(n, 3)))
op = CoulombOperator([bg_i], rng)
E0 = ens.kinetic_energy_eV().mean()
tau = slowing_down_time(DEUTERON, 1.5e3, bg_i)
dt = tau / 200
for _ in range(100):
    op(ens, np.arange(n), dt, 0.0)
out["collisions"]["maxwellian_energy_drift_over_half_tau"] = float(ens.kinetic_energy_eV().mean() / E0 - 1)

# ---------------------------------------------------------------- fusion
out["fusion"] = fusion_verify()
# beam-target MC vs n sigma v for 100 keV deuterons on cold tritons
rng = np.random.default_rng(2)
n = 5000
E_keV = 100.0
vD = math.sqrt(2 * E_keV * 1e3 * E_CHARGE / DEUTERON.mass_kg)
ens = Ensemble(DEUTERON, np.zeros((n, 3)), np.tile([vD, 0, 0], (n, 1)))
from ccsim.constants import TRITON
samp = FusionEventSampler("DT", TRITON, 1e26, 1.0, rng)   # dense cold target so the count is large
dt = 1e-9
tot = 0.0
for _ in range(200):
    c = samp(ens, np.arange(n), dt, 0.0)
    tot += c["fusion_events"]
mu = DEUTERON.mass_kg * TRITON.mass_kg / (DEUTERON.mass_kg + TRITON.mass_kg)
E_cm = 0.5 * mu * vD**2 / (1e3 * E_CHARGE)
expected = n * 200 * dt * 1e26 * float(cross_section_m2("DT", E_cm)[0]) * vD
out["fusion"]["beam_target_mc_vs_n_sigma_v_rel"] = float(tot / expected - 1)
out["fusion"]["beam_target_expected_events"] = expected

# ---------------------------------------------------------------- chemistry
w_ch4 = imag_freq_from_width(0.815, -0.02, 0.55, 0.94)
out["chemistry"] = {
    "eckart_ch4_imag_freq_cm1": w_ch4, "eckart_ch4_Tc_K": crossover_temperature_K(w_ch4),
    "eckart_gamma_900K": tunnelling_factor(900, 0.815, -0.02, w_ch4), "eckart_gamma_300K": tunnelling_factor(300, 0.815, -0.02, w_ch4),
    "qsurface_reference": {"imag_freq_cm1": 1256.1, "Tc_K": 287.6, "gamma_900K": 1.2030, "gamma_300K": 7.662},
}
# EEDF integration: constant sigma -> k = sigma * <v> = sigma * sqrt(8 kTe / (pi m_e))
Te = 3.0
k_num = eedf_rate_coefficient(lambda E: np.full_like(E, 1e-20), Te, E_max_eV=200.0, n=20000)
k_th = 1e-20 * math.sqrt(8 * Te * E_CHARGE / (math.pi * M_E))
out["chemistry"]["eedf_constant_sigma_rel_err"] = float(abs(k_num / k_th - 1))
net = co2_h2_plasma_network()
sn = sabatier_surface_network(0.05, 0.3, 0.1, 600.0)
out["chemistry"]["networks_balanced"] = {"plasma": len(net.reactions), "surface": len(sn.reactions)}

# ---- sphere winding: turns equally spaced in z give K ∝ sin θ → exactly uniform interior field ----
from ccsim.geometry import sphere_winding, close_with_return, Winding
from ccsim.fields import biot_savart
from ccsim.constants import MU0
_sw, _na = close_with_return(sphere_winding(0.76, 24), 3.8)
_w = Winding("sphere", _sw, 1.0, 0.005, closed=True, active_count=_na)
_P = np.array([[0, 0, 0], [0.25, 0.1, -0.15], [-0.3, 0.2, 0.3], [0.1, -0.35, 0.05], [0.0, 0.0, 0.45]])
_B = biot_savart(_P, _w)
_c = 0.94  # the winding stops at |z| = 0.94 R: K = K0 sin θ on θ ∈ [θ0, π−θ0], K0 = N I / (2 c R), B(0) = μ0 K0 (c − c³/3)
_B_th = MU0 * 24 * 1.0 / (2 * _c * 0.76) * (_c - _c**3 / 3)
out["fields_sphere_winding"] = {"B_centre_rel_err": float(abs(_B[0, 2] / _B_th - 1)), "B_axial_rel_err_max": float(np.max(np.abs(_B[:, 2] / _B_th - 1))),
                                "B_transverse_over_axial_max": float(np.max(np.hypot(_B[:, 0], _B[:, 1]) / _B[:, 2])),
                                "B_theory_T_per_A": _B_th}

# ---- fast kernel (numba or numpy fallback) vs the reference Biot–Savart ----
from ccsim.fastfield import SegmentSource, HAVE_NUMBA
_src = SegmentSource(_w)
_Pk = np.random.default_rng(3).uniform(-0.6, 0.6, (40, 3))
_Bk, _Br = _src.B(_Pk), biot_savart(_Pk, _w)
out["fastfield"] = {"numba": HAVE_NUMBA, "rel_err_max": float(np.abs(_Bk - _Br).max() / np.abs(_Br).max())}

# ---- numba orbit integrator vs run_orbits (same particles, same field) ----
try:
    from ccsim.fastorbits import HAVE_NUMBA as _HN, fast_orbits as _fo
    from ccsim.geometry import codex_toroidal as _ct, assembly_windings as _aw, MultiWinding as _MW
    from ccsim.fields import FieldGrid as _FG
    from ccsim.particles import monoenergetic_ensemble as _me, run_orbits as _ro, sample_positions_in_vessel as _sp
    from ccsim.presets import PRESETS as _PR
    from ccsim.runner import make_field as _mf, roi_mask as _rm
    from ccsim.constants import DEUTERON as _D, E_CHARGE as _E
    _pr = _PR["reactor"](); _w = _MW("p", _aw(_ct("precess", 0.22, 24.0, 3.0), "pair-side-crossed", 1.0, 0.005, 1.0, "precess", core_scale=0.45))
    _g = _FG(_w, 0.87, 17); _pool = _sp(400, 0.18, 0.58, np.random.default_rng(271828)); _vD = math.sqrt(2 * 15e3 * _E / _D.mass_kg)
    _f, _a = _mf(_g, _pr, _w, _pr.current_A, _pr.wire_radius_m); _brms = float(abs(_f.gain) * np.sqrt(np.mean(_g.Bmag.ravel()[_rm(_g)] ** 2)))
    _dt = min(_pr.dt_s, 0.15 / (_E / _D.mass_kg * _brms)); _pos = _pool[_g.wire_distance_at(_pool) > 0.03][:48]
    _e1 = _me(_D, 48, 15e3, _pos, np.random.default_rng(1)); _e2 = _me(_D, 48, 15e3, _pos, np.random.default_rng(1))
    _r = _ro(_e1, _f, 6 * 0.82 / _vD, _dt, 0.82); _lt, _k = _fo(_e2, _f, 6 * 0.82 / _vD, _dt, 0.82)
    out["fastorbits"] = {"numba": _HN, "retained_run_orbits": int(_r.retained), "retained_fast": int(np.sum(~np.isfinite(_lt))),
                         "S4_run_orbits": float(1 - np.sum(_r.loss_times <= 4 * 0.82 / _vD) / 48), "S4_fast": float(np.mean(~np.isfinite(_lt) | (_lt > 4 * 0.82 / _vD)))}
except Exception as _ex:  # pragma: no cover
    out["fastorbits"] = {"error": str(_ex)}

out["elapsed_s"] = time.time() - t_start
(ROOT / "results").mkdir(exist_ok=True)
with open(ROOT / "results" / "exp01_verify.json", "w") as fh:
    json.dump(out, fh, indent=2, default=float)
print(json.dumps(out, indent=1, default=float))
