"""Reaction networks for the gas, the plasma electrons and the wall.

Components
----------
* `Reaction` / `Network`: species with phases (gas, electron, surface, site),
  stoichiometric reactions with a rate law:
    - arrhenius        k = A · T^n · exp(−Ea/kT) · Γ_tunnel(T)   (Γ from an Eckart
                       barrier given V0, dE, ω‡ — the three numbers a CI-NEB gives)
    - electron_impact  k(Te) = ∫ σ(E) v f_Maxwell(E; Te) dE   (σ from a model or a table)
    - constant         k given
    - flux             surface arrival flux × sticking probability (from the
                       particle-tracing wall-flux map)
* Solvers: deterministic mass-action ODE (stiff, Radau) on number densities,
  and a Gillespie SSA for small populations (surface patches).
* Electron-impact cross-section models: Lotz ionisation (semi-empirical,
  decent to ~30 %), threshold-linear dissociation (**placeholder**: shape only),
  and tabulated (E, σ) arrays for LXCat-style data.
* `load_qsurface_network` imports a Q-Surface reaction-network JSON and turns
  every sech2 barrier into (V0, dE, ω‡) with exact Eckart tunnelling.
* `WallFluxMap` bins particle-tracing losses into ion flux per wall patch, the
  quantity that drives a surface network's collision step.

Provenance rule: every preset reaction carries a `source` string.  Values
marked PLACEHOLDER are shape-only and must be replaced before quantitative use.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import json
import math
from typing import Callable, Dict, List, Mapping, Optional, Sequence, Tuple

import numpy as np
from scipy.integrate import solve_ivp

from .constants import E_CHARGE, K_B, M_E
from .eckart import imag_freq_from_width, tunnelling_factor

KB_EV = K_B / E_CHARGE


# ---------------------------------------------------------------------------
# Electron-impact cross-sections and EEDF-averaged rate coefficients
# ---------------------------------------------------------------------------


def lotz_ionization_m2(E_eV, ionization_eV: float, electrons_in_shell: int = 2, a: float = 4.5e-18) -> np.ndarray:
    """Lotz (1968) single-shell ionisation cross-section: σ = a·q·ln(E/I)/(E·I)  [m²], a in eV²·m²."""

    E = np.atleast_1d(np.asarray(E_eV, dtype=float))
    with np.errstate(divide="ignore", invalid="ignore"):
        s = a * electrons_in_shell * np.log(np.maximum(E / ionization_eV, 1.0)) / (E * ionization_eV)
    return np.where(E > ionization_eV, s, 0.0)   # a = 4.5e-14 cm² eV² = 4.5e-18 m² eV² (Lotz's universal constant)


def threshold_linear_m2(E_eV, threshold_eV: float, sigma_max_m2: float, E_peak_eV: float) -> np.ndarray:
    """PLACEHOLDER shape: rises linearly from threshold to σ_max at E_peak, then falls as 1/E."""

    E = np.atleast_1d(np.asarray(E_eV, dtype=float))
    rise = sigma_max_m2 * (E - threshold_eV) / max(E_peak_eV - threshold_eV, 1e-9)
    fall = sigma_max_m2 * E_peak_eV / np.maximum(E, 1e-9)
    return np.where(E <= threshold_eV, 0.0, np.where(E <= E_peak_eV, rise, fall))


def tabulated_m2(table: Sequence[Tuple[float, float]]) -> Callable[[np.ndarray], np.ndarray]:
    """Cross-section from an (E [eV], σ [m²]) table (e.g. LXCat export), log-linear interpolation."""

    E_t = np.array([row[0] for row in table], dtype=float)
    s_t = np.array([row[1] for row in table], dtype=float)

    def f(E_eV):
        E = np.atleast_1d(np.asarray(E_eV, dtype=float))
        return np.where(E < E_t[0], 0.0, np.interp(E, E_t, s_t, right=s_t[-1] * E_t[-1] / np.maximum(E, 1e-9)))

    return f


def eedf_rate_coefficient(sigma: Callable[[np.ndarray], np.ndarray], Te_eV: float, E_max_eV: float = 300.0, n: int = 4000) -> float:
    """k(Te) = ∫ σ(E) v(E) f_M(E) dE with a Maxwellian EEDF (m³/s)."""

    E = np.linspace(0.0, E_max_eV, n)
    v = np.sqrt(2.0 * E * E_CHARGE / M_E)
    f = 2.0 * np.sqrt(E / math.pi) * Te_eV ** (-1.5) * np.exp(-E / Te_eV)
    return float(np.trapezoid(sigma(E) * v * f, E))


# ---------------------------------------------------------------------------
# Network definition
# ---------------------------------------------------------------------------


@dataclass
class SpeciesDef:
    name: str
    phase: str = "gas"           # gas | electron | surface | site
    composition: Dict[str, int] = field(default_factory=dict)
    sites: int = 0


@dataclass
class Reaction:
    id: str
    reactants: Dict[str, int]
    products: Dict[str, int]
    law: str                                   # arrhenius | electron_impact | constant | flux
    A: float = 0.0                             # prefactor (units depend on order)
    n: float = 0.0
    Ea_eV: float = 0.0
    tunnel: Optional[Tuple[float, float, float]] = None   # (V0, dE, imag_freq_cm1)
    sigma: Optional[Callable] = None           # electron_impact cross-section E[eV] -> m^2
    k_const: float = 0.0
    flux: float = 0.0                          # for 'flux' law: arrival rate coefficient
    sticking: float = 1.0
    source: str = "unspecified"
    note: str = ""

    def rate_coefficient(self, T_K: float, Te_eV: Optional[float] = None) -> float:
        if self.law == "arrhenius":
            k = self.A * T_K**self.n * math.exp(-self.Ea_eV / (KB_EV * T_K))
            if self.tunnel is not None:
                k *= tunnelling_factor(T_K, *self.tunnel)
            return k
        if self.law == "electron_impact":
            if Te_eV is None:
                raise ValueError(f"{self.id}: electron temperature needed")
            return eedf_rate_coefficient(self.sigma, Te_eV)
        if self.law == "constant":
            return self.k_const
        if self.law == "flux":
            return self.flux * self.sticking
        raise ValueError(f"unknown rate law {self.law}")


@dataclass
class Network:
    name: str
    species: List[SpeciesDef]
    reactions: List[Reaction]

    def __post_init__(self):
        names = [s.name for s in self.species]
        if len(names) != len(set(names)):
            raise ValueError("duplicate species")
        self.index = {n: i for i, n in enumerate(names)}
        for r in self.reactions:
            for side in (r.reactants, r.products):
                for n in side:
                    if n not in self.index:
                        raise ValueError(f"{r.id}: unknown species {n}")
            self._check_balance(r)

    def _check_balance(self, r: Reaction):
        comp = {s.name: s.composition for s in self.species}
        sites = {s.name: s.sites for s in self.species}
        lhs, rhs = {}, {}
        for n, c in r.reactants.items():
            for el, k in comp[n].items():
                lhs[el] = lhs.get(el, 0) + c * k
        for n, c in r.products.items():
            for el, k in comp[n].items():
                rhs[el] = rhs.get(el, 0) + c * k
        if lhs != rhs:
            raise ValueError(f"{r.id}: element balance fails {lhs} -> {rhs}")
        ls = sum(c * sites[n] for n, c in r.reactants.items())
        rs = sum(c * sites[n] for n, c in r.products.items())
        if ls != rs:
            raise ValueError(f"{r.id}: site balance fails {ls} -> {rs}")

    def _kinetic_orders(self):
        """Reactant (index, order) pairs that enter the rate; gas reactants of 'flux' steps are reservoirs."""

        phase = {s.name: s.phase for s in self.species}
        out = []
        for r in self.reactions:
            terms = [(self.index[n], c) for n, c in r.reactants.items() if not (r.law == "flux" and phase[n] == "gas")]
            out.append(terms)
        return out

    def stoichiometry(self) -> np.ndarray:
        S = np.zeros((len(self.species), len(self.reactions)))
        for j, r in enumerate(self.reactions):
            for n, c in r.reactants.items():
                S[self.index[n], j] -= c
            for n, c in r.products.items():
                S[self.index[n], j] += c
        return S

    def rate_coefficients(self, T_K: float, Te_eV: Optional[float] = None) -> np.ndarray:
        return np.array([r.rate_coefficient(T_K, Te_eV) for r in self.reactions])

    # -- deterministic mass action ------------------------------------------
    def integrate(self, n0: Mapping[str, float], t_end: float, T_K: float, Te_eV: Optional[float] = None,
                  n_out: int = 201, fixed: Sequence[str] = (), rtol: float = 1e-7, atol: float = 1e-12):
        """Mass-action ODE on number densities (m^-3); `fixed` species are held constant (reservoirs)."""

        S = self.stoichiometry()
        k = self.rate_coefficients(T_K, Te_eV)
        y0 = np.zeros(len(self.species))
        for name, val in n0.items():
            y0[self.index[name]] = val
        fixed_idx = [self.index[f] for f in fixed]
        orders = self._kinetic_orders()

        def rhs(t, y):
            y = np.maximum(y, 0.0)
            rates = k.copy()
            for j, terms in enumerate(orders):
                for i, c in terms:
                    rates[j] *= y[i] ** c
            dy = S @ rates
            dy[fixed_idx] = 0.0
            return dy

        t_eval = np.linspace(0.0, t_end, n_out)
        sol = solve_ivp(rhs, (0.0, t_end), y0, t_eval=t_eval, method="Radau", rtol=rtol, atol=atol)
        return sol.t, {s.name: sol.y[i] for i, s in enumerate(self.species)}, k

    # -- Gillespie on integer populations ------------------------------------
    def gillespie(self, counts0: Mapping[str, int], t_end: float, T_K: float, rng, volume_m3: float = 1.0,
                  Te_eV: Optional[float] = None, max_events: int = 2_000_000, n_out: int = 201, fixed: Sequence[str] = ()):
        """Stochastic simulation on molecule counts; k for bimolecular steps is divided by the volume."""

        S = self.stoichiometry()
        k = self.rate_coefficients(T_K, Te_eV)
        x = np.zeros(len(self.species), dtype=np.int64)
        for name, val in counts0.items():
            x[self.index[name]] = int(val)
        fixed_idx = set(self.index[f] for f in fixed)
        orders = self._kinetic_orders()
        order_total = [sum(c for _, c in terms) for terms in orders]
        t = 0.0
        t_grid = np.linspace(0.0, t_end, n_out)
        hist = np.zeros((n_out, len(self.species)))
        gi = 0
        events = 0
        while events < max_events:
            while gi < n_out and t_grid[gi] <= t:
                hist[gi] = x
                gi += 1
            a = np.empty(len(self.reactions))
            for j, terms in enumerate(orders):
                prop = k[j] / volume_m3 ** max(order_total[j] - 1, 0)
                for i, c in terms:
                    for m in range(c):
                        prop *= max(x[i] - m, 0)
                    if c > 1:
                        prop /= math.factorial(c)
                a[j] = prop
            total = a.sum()
            if total <= 0:
                break
            t += -math.log(max(rng.random(), 1e-300)) / total
            if t > t_end:
                break
            j = int(np.searchsorted(np.cumsum(a), rng.random() * total, side="right"))
            j = min(j, len(a) - 1)
            for i in range(len(x)):
                if i in fixed_idx:
                    continue
                x[i] += int(S[i, j])
            events += 1
        while gi < n_out:
            hist[gi] = x
            gi += 1
        return t_grid, {s.name: hist[:, i] for i, s in enumerate(self.species)}, events


# ---------------------------------------------------------------------------
# Wall-flux coupling from particle tracing
# ---------------------------------------------------------------------------


@dataclass
class WallFluxMap:
    """Ion arrival flux on a spherical wall binned in (cosθ, φ) patches."""

    n_theta: int
    n_phi: int
    flux_m2_s: np.ndarray       # (n_theta, n_phi)
    total_rate_s: float
    wall_radius_m: float

    @classmethod
    def from_losses(cls, loss_positions: np.ndarray, loss_kinds: np.ndarray, sim_time_s: float, weight: float,
                    wall_radius_m: float, n_theta: int = 12, n_phi: int = 24) -> "WallFluxMap":
        P = np.asarray(loss_positions)
        mask = np.asarray(loss_kinds) == "wall"
        P = P[mask]
        flux = np.zeros((n_theta, n_phi))
        if len(P):
            r = np.linalg.norm(P, axis=1)
            ct = P[:, 2] / np.maximum(r, 1e-300)
            ph = np.arctan2(P[:, 1], P[:, 0])
            it = np.clip(((ct + 1.0) / 2.0 * n_theta).astype(int), 0, n_theta - 1)
            ip = np.clip(((ph + math.pi) / (2 * math.pi) * n_phi).astype(int), 0, n_phi - 1)
            np.add.at(flux, (it, ip), 1.0)
        patch_area = 4.0 * math.pi * wall_radius_m**2 / (n_theta * n_phi)   # equal-area in cosθ
        flux = flux * weight / (sim_time_s * patch_area)
        return cls(n_theta, n_phi, flux, float(len(P) * weight / sim_time_s), wall_radius_m)

    def peak_to_mean(self) -> float:
        m = self.flux_m2_s.mean()
        return float(self.flux_m2_s.max() / m) if m > 0 else float("nan")

    def monolayers_per_s(self, site_density_m2: float = 1.5e19) -> np.ndarray:
        """Arrival flux in ML/s per patch (Pt(111) has 1.5e19 sites/m²)."""

        return self.flux_m2_s / site_density_m2


# ---------------------------------------------------------------------------
# Q-Surface import
# ---------------------------------------------------------------------------


def load_qsurface_network(path: str, flux_ml_per_s: float = 0.0, sticking: float = 0.0) -> Network:
    """Import a Q-Surface network JSON (as written by `qsurf export-preset`).

    Thermal steps become Arrhenius(prefactor·degeneracy·steric, Ea = height) with an
    exact Eckart tunnelling factor built from (height, reaction_energy, ω‡) where
    ω‡ is *derived from the sech2 shape* (only m·w² matters).  The collision step
    becomes a 'flux' law whose sticking must be supplied (from the beam / thermal
    gas average) — or from the wall-flux map of a particle run.
    """

    with open(path) as fh:
        data = json.load(fh)
    species = [SpeciesDef(s["name"], s["phase"], dict(s.get("composition", {})), int(s.get("sites", 0))) for s in data["species"]]
    species.append(SpeciesDef("*", "site", {}, 1))
    rxns = []
    for st in data["steps"]:
        b = st["barrier"]
        V0, dE = float(b["height_eV"]), float(b.get("reaction_energy_eV", 0.0))
        if b.get("model", "sech2") == "sech2":
            w = imag_freq_from_width(V0, dE, float(b["width_angstrom"]), float(b["effective_mass_amu"]))
            tunnel = (V0, dE, w)
        elif b.get("model") == "parabolic":
            tunnel = (V0, dE, float(b.get("imaginary_frequency_cm1", 0.0))) if b.get("imaginary_frequency_cm1", 0) > 0 else None
        else:
            tunnel = None
        reactants = {k: int(v) for k, v in st["reactants"].items()}
        products = {k: int(v) for k, v in st["products"].items()}
        if st["kind"] == "thermal":
            rxns.append(Reaction(st["id"], reactants, products, "arrhenius",
                                 A=float(st.get("prefactor_s", 1e13)) * float(st.get("degeneracy", 1.0)) * float(st.get("steric_factor", 1.0)),
                                 Ea_eV=V0, tunnel=tunnel, source=st.get("provenance", "qsurf"), note=st.get("calibration_status", "")))
        else:
            gas = [n for n in reactants if n != "*" and any(s.name == n and s.phase == "gas" for s in species)]
            rxns.append(Reaction(st["id"], reactants, products, "flux",
                                 flux=flux_ml_per_s, sticking=sticking, tunnel=tunnel, Ea_eV=V0,
                                 source=st.get("provenance", "qsurf"), note="gas reactant treated as reservoir; flux × sticking"))
    return Network(data.get("name", "qsurface"), species, rxns)


# ---------------------------------------------------------------------------
# Presets
# ---------------------------------------------------------------------------


def co2_h2_plasma_network(Te_eV: float = 3.0) -> Network:
    """Reduced CO₂/H₂ plasma-chemistry network (gas phase + electron impact).

    Sources: O + H2 and OH + H2 rate expressions are the standard combustion
    values (Baulch et al. evaluation, T^n form, k in m³/s); CO + O + M and
    H + H + M three-body rates are order-of-magnitude literature values;
    electron-impact dissociation cross-sections are PLACEHOLDER shapes with
    thresholds at the known excitation/dissociation energies; ionisation uses
    Lotz.  This network is a scaffold for LXCat data, not a validated model.
    """

    species = [SpeciesDef("e", "electron"), SpeciesDef("CO2", "gas", {"C": 1, "O": 2}), SpeciesDef("CO", "gas", {"C": 1, "O": 1}),
               SpeciesDef("O", "gas", {"O": 1}), SpeciesDef("O2", "gas", {"O": 2}), SpeciesDef("H2", "gas", {"H": 2}),
               SpeciesDef("H", "gas", {"H": 1}), SpeciesDef("OH", "gas", {"O": 1, "H": 1}), SpeciesDef("H2O", "gas", {"H": 2, "O": 1}),
               SpeciesDef("CO2+", "gas", {"C": 1, "O": 2}), SpeciesDef("H2+", "gas", {"H": 2}), SpeciesDef("M", "gas")]
    rx = [
        Reaction("e_CO2_diss", {"e": 1, "CO2": 1}, {"e": 1, "CO": 1, "O": 1}, "electron_impact",
                 sigma=lambda E: threshold_linear_m2(E, 7.0, 1.0e-20, 12.0), source="PLACEHOLDER shape; threshold ~7 eV (CO2 electronic dissociation)"),
        Reaction("e_H2_diss", {"e": 1, "H2": 1}, {"e": 1, "H": 2}, "electron_impact",
                 sigma=lambda E: threshold_linear_m2(E, 8.8, 6.0e-21, 15.0), source="PLACEHOLDER shape; threshold 8.8 eV (b³Σ dissociation)"),
        Reaction("e_CO2_ion", {"e": 1, "CO2": 1}, {"e": 2, "CO2+": 1}, "electron_impact",
                 sigma=lambda E: lotz_ionization_m2(E, 13.8, 4), source="Lotz formula, I = 13.8 eV"),
        Reaction("e_H2_ion", {"e": 1, "H2": 1}, {"e": 2, "H2+": 1}, "electron_impact",
                 sigma=lambda E: lotz_ionization_m2(E, 15.4, 2), source="Lotz formula, I = 15.4 eV"),
        Reaction("O_H2", {"O": 1, "H2": 1}, {"OH": 1, "H": 1}, "arrhenius", A=3.44e-19 / 300.0**2.67, n=2.67, Ea_eV=3160.0 * KB_EV,
                 source="Baulch et al.: k = 3.44e-13 (T/300)^2.67 exp(-3160/T) cm³/s"),
        Reaction("OH_H2", {"OH": 1, "H2": 1}, {"H2O": 1, "H": 1}, "arrhenius", A=7.7e-18, Ea_eV=2100.0 * KB_EV,
                 source="Baulch et al.: k = 7.7e-12 exp(-2100/T) cm³/s"),
        Reaction("CO_O_M", {"CO": 1, "O": 1, "M": 1}, {"CO2": 1, "M": 1}, "arrhenius", A=1.7e-45, Ea_eV=1510.0 * KB_EV,
                 source="order-of-magnitude three-body value (~1.7e-33 exp(-1510/T) cm⁶/s)"),
        Reaction("H_H_M", {"H": 2, "M": 1}, {"H2": 1, "M": 1}, "arrhenius", A=6.0e-45, n=-1.0, Ea_eV=0.0,
                 source="order-of-magnitude three-body value (~6e-33 (T/300)^-1 cm⁶/s)"),
        Reaction("O_O_M", {"O": 2, "M": 1}, {"O2": 1, "M": 1}, "arrhenius", A=5.0e-46, Ea_eV=0.0,
                 source="order-of-magnitude three-body value (~5e-34 cm⁶/s)"),
    ]
    return Network("CO2/H2 reduced plasma network", species, rx)


def sabatier_surface_network(flux_ml_per_s: float, sticking_co2: float, sticking_h2: float, T_K: float) -> Network:
    """Illustrative surface methanation network on a catalyst patch (site balance enforced).

    CO2* + H* → COOH*, COOH* → CO* + OH*, CO* + H* → CHO*, ... down to CH4(g) and
    H2O(g).  Barriers are ILLUSTRATIVE (0.6–1.1 eV, typical of DFT ranges on
    Ni/Pt) and are the fields your MACE chain is meant to replace: each step's
    (V0, dE, ω‡) plugs straight into `tunnel`.
    """

    S = SpeciesDef
    species = [S("CO2(g)", "gas", {"C": 1, "O": 2}), S("H2(g)", "gas", {"H": 2}), S("CH4(g)", "gas", {"C": 1, "H": 4}), S("H2O(g)", "gas", {"H": 2, "O": 1}),
               S("CO2*", "surface", {"C": 1, "O": 2}, 1), S("H*", "surface", {"H": 1}, 1), S("COOH*", "surface", {"C": 1, "O": 2, "H": 1}, 1),
               S("CO*", "surface", {"C": 1, "O": 1}, 1), S("OH*", "surface", {"O": 1, "H": 1}, 1), S("C*", "surface", {"C": 1}, 1),
               S("O*", "surface", {"O": 1}, 1), S("CH*", "surface", {"C": 1, "H": 1}, 1), S("CH2*", "surface", {"C": 1, "H": 2}, 1),
               S("CH3*", "surface", {"C": 1, "H": 3}, 1), S("*", "site", {}, 1)]

    def ar(id_, r, p, V0, dE, w, A=1e13):
        return Reaction(id_, r, p, "arrhenius", A=A, Ea_eV=V0, tunnel=(V0, dE, w), source="ILLUSTRATIVE barrier; replace with CI-NEB (V0, dE, ω‡)")

    rx = [
        Reaction("CO2_ads", {"CO2(g)": 1, "*": 1}, {"CO2*": 1}, "flux", flux=flux_ml_per_s, sticking=sticking_co2, source="wall-flux × sticking"),
        Reaction("H2_ads", {"H2(g)": 1, "*": 2}, {"H*": 2}, "flux", flux=flux_ml_per_s, sticking=sticking_h2, source="wall-flux × sticking (dissociative)"),
        ar("COOH_form", {"CO2*": 1, "H*": 1}, {"COOH*": 1, "*": 1}, 0.75, 0.10, 1200.0),
        ar("COOH_split", {"COOH*": 1, "*": 1}, {"CO*": 1, "OH*": 1}, 0.60, -0.30, 900.0),
        ar("CO_diss", {"CO*": 1, "*": 1}, {"C*": 1, "O*": 1}, 1.10, 0.30, 800.0),
        ar("C_H", {"C*": 1, "H*": 1}, {"CH*": 1, "*": 1}, 0.70, -0.20, 1400.0),
        ar("CH_H", {"CH*": 1, "H*": 1}, {"CH2*": 1, "*": 1}, 0.65, 0.05, 1350.0),
        ar("CH2_H", {"CH2*": 1, "H*": 1}, {"CH3*": 1, "*": 1}, 0.60, -0.05, 1300.0),
        ar("CH3_H", {"CH3*": 1, "H*": 1}, {"CH4(g)": 1, "*": 2}, 0.80, -0.10, 1250.0),
        ar("O_H", {"O*": 1, "H*": 1}, {"OH*": 1, "*": 1}, 0.90, 0.10, 1000.0),
        ar("OH_H", {"OH*": 1, "H*": 1}, {"H2O(g)": 1, "*": 2}, 0.85, -0.20, 1100.0),
    ]
    return Network("illustrative Sabatier surface network", species, rx)
