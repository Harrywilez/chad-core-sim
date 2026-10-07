"""numba fast path for the vacuum-field Boris integration used by the design evaluator.

Same physics as `particles.run_orbits` for the collisionless, E = 0 case: relativistic
Boris, trilinear B from the unit-scale grid with gain I/(I_grid·scale), per-particle
substepping so that Ω·dt_sub ≤ omega_dt_max, spherical wall, loss on leaving the grid or
coming within the conductor radius (checked against the grid's conductor-distance map,
refined with the exact strand distance when close).  Returns per-particle loss times
(inf = retained) and loss kinds.  Verified against run_orbits in exp01 (`fastorbits`).
"""

from __future__ import annotations

import math

import numpy as np

from .constants import C_LIGHT

try:
    from numba import njit
    HAVE_NUMBA = True
except Exception:  # pragma: no cover
    HAVE_NUMBA = False

    def njit(*a, **k):
        def deco(f):
            return f
        return deco if not (a and callable(a[0])) else a[0]


@njit(cache=True, fastmath=True)
def _trilinear(G, D, half, dx, n, x, y, z, out):
    """B (3) and conductor distance at (x,y,z) from the grid; returns distance."""
    fx = (x + half) / dx; fy = (y + half) / dx; fz = (z + half) / dx
    ix = int(math.floor(fx)); iy = int(math.floor(fy)); iz = int(math.floor(fz))
    if ix < 0: ix = 0
    if iy < 0: iy = 0
    if iz < 0: iz = 0
    if ix > n - 2: ix = n - 2
    if iy > n - 2: iy = n - 2
    if iz > n - 2: iz = n - 2
    wx = fx - ix; wy = fy - iy; wz = fz - iz
    if wx < 0: wx = 0.0
    if wy < 0: wy = 0.0
    if wz < 0: wz = 0.0
    if wx > 1: wx = 1.0
    if wy > 1: wy = 1.0
    if wz > 1: wz = 1.0
    out[0] = 0.0; out[1] = 0.0; out[2] = 0.0
    d = 0.0
    for a in range(2):
        ax = wx if a else 1.0 - wx
        for b in range(2):
            ay = wy if b else 1.0 - wy
            for c in range(2):
                az = wz if c else 1.0 - wz
                w = ax * ay * az
                out[0] += w * G[ix + a, iy + b, iz + c, 0]
                out[1] += w * G[ix + a, iy + b, iz + c, 1]
                out[2] += w * G[ix + a, iy + b, iz + c, 2]
                d += w * D[ix + a, iy + b, iz + c]
    return d


@njit(cache=True, fastmath=True)
def _fast_orbits(x, v, q_over_m, G, D, half, dx, n, gain, scale, dt, t_max, wall, omega_dt_max, wire_radius, wire_pts, near_thresh, max_sub):
    N = x.shape[0]
    loss_t = np.full(N, np.inf)
    kind = np.zeros(N, np.int8)   # 0 retained, 1 wall, 2 wire, 3 grid-exit, 4 numerical
    steps = int(math.ceil(t_max / dt))
    c2 = C_LIGHT * C_LIGHT
    B = np.zeros(3)
    half_m = half * scale
    for p in range(N):
        px, py, pz = x[p, 0], x[p, 1], x[p, 2]
        vx, vy, vz = v[p, 0], v[p, 1], v[p, 2]
        t = 0.0
        for s in range(steps):
            d = _trilinear(G, D, half, dx, n, px / scale, py / scale, pz / scale, B)
            bx, by, bz = gain * B[0], gain * B[1], gain * B[2]
            bn = math.sqrt(bx * bx + by * by + bz * bz)
            nsub = int(math.ceil(abs(q_over_m) * bn * dt / omega_dt_max))
            if nsub < 1: nsub = 1
            if nsub > max_sub: nsub = max_sub
            sdt = dt / nsub
            for k in range(nsub):
                if k > 0:
                    d = _trilinear(G, D, half, dx, n, px / scale, py / scale, pz / scale, B)
                    bx, by, bz = gain * B[0], gain * B[1], gain * B[2]
                # relativistic Boris, E = 0
                v2 = vx * vx + vy * vy + vz * vz
                r = v2 / c2
                if r > 0.999999: r = 0.999999
                gam = 1.0 / math.sqrt(1.0 - r)
                ux, uy, uz = gam * vx, gam * vy, gam * vz
                gm = math.sqrt(1.0 + (ux * ux + uy * uy + uz * uz) / c2)
                f = 0.5 * q_over_m * sdt / gm
                tx, ty, tz = f * bx, f * by, f * bz
                t2 = tx * tx + ty * ty + tz * tz
                sx, sy, sz = 2.0 * tx / (1.0 + t2), 2.0 * ty / (1.0 + t2), 2.0 * tz / (1.0 + t2)
                upx = ux + (uy * tz - uz * ty); upy = uy + (uz * tx - ux * tz); upz = uz + (ux * ty - uy * tx)
                ux2 = ux + (upy * sz - upz * sy); uy2 = uy + (upz * sx - upx * sz); uz2 = uz + (upx * sy - upy * sx)
                gn = math.sqrt(1.0 + (ux2 * ux2 + uy2 * uy2 + uz2 * uz2) / c2)
                vx, vy, vz = ux2 / gn, uy2 / gn, uz2 / gn
                px += vx * sdt; py += vy * sdt; pz += vz * sdt
            t += dt
            rr = math.sqrt(px * px + py * py + pz * pz)
            if not (rr == rr):
                loss_t[p] = t; kind[p] = 4; break
            if rr >= wall:
                loss_t[p] = t; kind[p] = 1; break
            if abs(px) >= half_m or abs(py) >= half_m or abs(pz) >= half_m:
                loss_t[p] = t; kind[p] = 3; break
            # conductor proximity: interpolated distance map (unit scale), refine with the exact strand points
            d = _trilinear(G, D, half, dx, n, px / scale, py / scale, pz / scale, B) * scale
            if d <= near_thresh:
                best = 1e300
                for q in range(wire_pts.shape[0]):
                    ddx = wire_pts[q, 0] - px; ddy = wire_pts[q, 1] - py; ddz = wire_pts[q, 2] - pz
                    dd = ddx * ddx + ddy * ddy + ddz * ddz
                    if dd < best: best = dd
                if math.sqrt(best) <= wire_radius:
                    loss_t[p] = t; kind[p] = 2; break
        x[p, 0], x[p, 1], x[p, 2] = px, py, pz
        v[p, 0], v[p, 1], v[p, 2] = vx, vy, vz
    return loss_t, kind


def fast_orbits(ens, field, t_max: float, dt: float, wall_radius: float, omega_dt_max: float = 0.2, max_substeps: int = 4096):
    """Drop-in for run_orbits in the collisionless E = 0 case. Returns (loss_times (N, inf=retained), kinds)."""
    g = field.grid
    q_over_m = ens.species.charge_C / ens.species.mass_kg
    x = np.ascontiguousarray(ens.x, dtype=float)
    v = np.ascontiguousarray(ens.v, dtype=float)
    near = max(3.0 * g.dx * field.scale, 2.0 * field.wire_radius)
    wire_pts = np.ascontiguousarray(field._wire_points, dtype=float)
    lt, kind = _fast_orbits(x, v, q_over_m, np.ascontiguousarray(g.B), np.ascontiguousarray(g.wire_distance), g.half, g.dx, g.n,
                            float(field.gain), float(field.scale), float(dt), float(t_max), float(wall_radius), float(omega_dt_max),
                            float(field.wire_radius), wire_pts, float(near), int(max_substeps))
    ens.x[:] = x; ens.v[:] = v
    names = np.array(["retained", "wall", "wire", "grid-exit", "numerical"])
    return lt, names[kind]
