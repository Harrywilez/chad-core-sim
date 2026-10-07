"""3-D file export for conductor geometries: centreline polylines (OBJ / CSV) and swept tube meshes
(binary STL, OBJ with one object per strand, binary PLY with per-strand vertex colours).

The tube sweep uses the rotation-minimising (Bishop) frame from `geometry.bishop_frame`; for a closed
strand the frame's holonomy is spread as a uniform twist along the curve so the last ring meets the
first one exactly and the mesh is watertight.  Open strands get flat end caps.

    from ccsim.export3d import tube_mesh, write_stl, write_obj_mesh, write_obj_polylines, write_ply, write_csv
"""

from __future__ import annotations

import math
import struct
from pathlib import Path
from typing import Iterable, Sequence, Tuple

import numpy as np

from .geometry import bishop_frame, unit


# ---------------------------------------------------------------------------
# Centrelines
# ---------------------------------------------------------------------------


def write_csv(path, points: np.ndarray, comment: str = "") -> None:
    """x,y,z rows (unit-ball units unless stated in the comment)."""
    p = np.asarray(points, dtype=float)
    with open(path, "w", encoding="utf-8") as fh:
        if comment:
            for line in comment.splitlines():
                fh.write(f"# {line}\n")
        fh.write("x,y,z\n")
        np.savetxt(fh, p, fmt="%.9f", delimiter=",")


def write_obj_polylines(path, curves: Sequence[Tuple[str, np.ndarray]], comment: str = "") -> None:
    """Polylines as OBJ `l` elements (one `o` object per curve, 2-vertex segments for maximum importer
    compatibility)."""
    with open(path, "w", encoding="utf-8") as fh:
        if comment:
            for line in comment.splitlines():
                fh.write(f"# {line}\n")
        base = 1
        for name, pts in curves:
            p = np.asarray(pts, dtype=float)
            fh.write(f"o {name}\n")
            fh.write("".join(f"v {x:.9f} {y:.9f} {z:.9f}\n" for x, y, z in p))
            fh.write("".join(f"l {base + i} {base + i + 1}\n" for i in range(len(p) - 1)))
            base += len(p)


# ---------------------------------------------------------------------------
# Tube meshes
# ---------------------------------------------------------------------------


def _periodic_frame(p: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Bishop frame on a closed polyline (p[0] == p[-1] not required; the last point connects to the
    first) with the holonomy removed by a uniform twist, so frame[n] == frame[0]."""
    n = len(p)
    ext = np.vstack((p[-2:], p, p[:3]))            # pad both ends so tangents are periodic
    t, n1, n2 = bishop_frame(ext)
    t, n1, n2 = t[2:2 + n + 1], n1[2:2 + n + 1], n2[2:2 + n + 1]    # frames at 0..n (n ≡ 0)
    # holonomy: angle from the transported normal at index n back to the normal at index 0 (same tangent)
    c = float(np.dot(n1[n], n1[0]))
    s = float(np.dot(np.cross(n1[n], n1[0]), t[0]))
    delta = math.atan2(s, c)
    ang = delta * np.arange(n + 1) / n
    ca, sa = np.cos(ang)[:, None], np.sin(ang)[:, None]
    m1 = ca * n1 + sa * n2
    m2 = -sa * n1 + ca * n2
    return t[:n], m1[:n], m2[:n]


def tube_mesh(points: np.ndarray, radius: float, n_seg: int = 8, closed: bool = False, caps: bool = True):
    """Sweep a circle of `radius` along the polyline.  Returns (V (M,3), F (K,3) int) — triangles with
    outward normals.  closed=True joins the last ring to the first (the duplicated closing point, if
    present, is dropped)."""
    p = np.asarray(points, dtype=float)
    if closed and np.linalg.norm(p[-1] - p[0]) < 1e-12:
        p = p[:-1]
    n = len(p)
    if closed:
        t, n1, n2 = _periodic_frame(p)
    else:
        t, n1, n2 = bishop_frame(p)
    ang = 2 * np.pi * np.arange(n_seg) / n_seg
    ca, sa = np.cos(ang), np.sin(ang)
    V = p[:, None, :] + radius * (ca[None, :, None] * n1[:, None, :] + sa[None, :, None] * n2[:, None, :])
    V = V.reshape(-1, 3)
    rings = n if closed else n - 1
    i = np.arange(rings)[:, None]
    j = np.arange(n_seg)[None, :]
    a = (i * n_seg + j).ravel()
    b = (i * n_seg + (j + 1) % n_seg).ravel()
    c = (((i + 1) % n) * n_seg + j).ravel()
    d = (((i + 1) % n) * n_seg + (j + 1) % n_seg).ravel()
    F = np.vstack((np.column_stack((a, c, b)), np.column_stack((b, c, d))))
    # orientation: make the normal of the first quad point away from the axis
    tri = V[F[0]]
    nrm = np.cross(tri[1] - tri[0], tri[2] - tri[0])
    if np.dot(nrm, V[F[0, 0]] - p[0]) < 0:
        F = F[:, [0, 2, 1]]
    if not closed and caps:
        base = len(V)
        V = np.vstack((V, p[[0]], p[[-1]]))
        ring0 = np.arange(n_seg)
        ringN = (n - 1) * n_seg + np.arange(n_seg)
        cap0 = np.column_stack((np.full(n_seg, base), ring0[(np.arange(n_seg) + 1) % n_seg], ring0))
        capN = np.column_stack((np.full(n_seg, base + 1), ringN, ringN[(np.arange(n_seg) + 1) % n_seg]))
        F = np.vstack((F, cap0, capN))
    return V, F.astype(np.int64)


def mesh_normals(V: np.ndarray, F: np.ndarray) -> np.ndarray:
    tri = V[F]
    n = np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0])
    return unit(n)


def write_stl(path, parts: Sequence[Tuple[str, np.ndarray, np.ndarray]], header: str = "") -> None:
    """Binary STL of several (name, V, F) parts (STL has no colours; all parts are merged)."""
    Fs = []
    Vs = []
    off = 0
    for _, V, F in parts:
        Vs.append(V)
        Fs.append(F + off)
        off += len(V)
    V = np.vstack(Vs)
    F = np.vstack(Fs)
    N = mesh_normals(V, F)
    tri = V[F]
    rec = np.zeros(len(F), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("attr", "<u2")])
    rec["n"] = N
    rec["v"] = tri
    with open(path, "wb") as fh:
        fh.write(header.encode("ascii", "replace")[:80].ljust(80, b" "))
        fh.write(struct.pack("<I", len(F)))
        fh.write(rec.tobytes())


def write_obj_mesh(path, parts: Sequence[Tuple[str, np.ndarray, np.ndarray]], comment: str = "") -> None:
    """OBJ mesh with one object (`o name`) per part and per-face normals omitted (viewers compute them)."""
    with open(path, "w", encoding="utf-8") as fh:
        if comment:
            for line in comment.splitlines():
                fh.write(f"# {line}\n")
        off = 1
        for name, V, F in parts:
            fh.write(f"o {name}\n")
            fh.write("".join(f"v {x:.6f} {y:.6f} {z:.6f}\n" for x, y, z in V))
            fh.write("".join(f"f {a + off} {b + off} {c + off}\n" for a, b, c in F))
            off += len(V)


def write_ply(path, parts: Sequence[Tuple[str, np.ndarray, np.ndarray]], colors: Sequence[Tuple[int, int, int]],
              comment: str = "") -> None:
    """Binary little-endian PLY with per-vertex RGB (one colour per part)."""
    Vs, Fs, Cs = [], [], []
    off = 0
    for (name, V, F), col in zip(parts, colors):
        Vs.append(V)
        Fs.append(F + off)
        Cs.append(np.tile(np.asarray(col, dtype=np.uint8), (len(V), 1)))
        off += len(V)
    V = np.vstack(Vs).astype("<f4")
    F = np.vstack(Fs).astype("<i4")
    C = np.vstack(Cs)
    vrec = np.zeros(len(V), dtype=[("x", "<f4"), ("y", "<f4"), ("z", "<f4"), ("r", "u1"), ("g", "u1"), ("b", "u1")])
    vrec["x"], vrec["y"], vrec["z"] = V[:, 0], V[:, 1], V[:, 2]
    vrec["r"], vrec["g"], vrec["b"] = C[:, 0], C[:, 1], C[:, 2]
    frec = np.zeros(len(F), dtype=[("n", "u1"), ("i", "<i4", 3)])
    frec["n"] = 3
    frec["i"] = F
    with open(path, "wb") as fh:
        head = ["ply", "format binary_little_endian 1.0"]
        head += [f"comment {line}" for line in comment.splitlines()]
        head += [f"element vertex {len(V)}", "property float x", "property float y", "property float z",
                 "property uchar red", "property uchar green", "property uchar blue",
                 f"element face {len(F)}", "property list uchar int vertex_indices", "end_header"]
        fh.write(("\n".join(head) + "\n").encode("ascii"))
        fh.write(vrec.tobytes())
        fh.write(frec.tobytes())


def mesh_is_closed(V: np.ndarray, F: np.ndarray) -> bool:
    """Every edge shared by exactly two faces (watertight, manifold)."""
    e = np.vstack((F[:, [0, 1]], F[:, [1, 2]], F[:, [2, 0]]))
    e = np.sort(e, axis=1)
    _, counts = np.unique(e, axis=0, return_counts=True)
    return bool(np.all(counts == 2))


def read_stl(path) -> Tuple[np.ndarray, np.ndarray]:
    """Binary STL → (normals, triangles (K,3,3))."""
    with open(path, "rb") as fh:
        fh.read(80)
        k = struct.unpack("<I", fh.read(4))[0]
        rec = np.frombuffer(fh.read(), dtype=[("n", "<f4", 3), ("v", "<f4", (3, 3)), ("attr", "<u2")], count=k)
    return rec["n"].copy(), rec["v"].copy()


def read_obj_polylines(path) -> dict:
    """OBJ `o`/`v`/`l` → {name: points (M,3)} (vertices in the order of first use along the `l` chain)."""
    verts = []
    out: dict = {}
    cur = None
    chain: list = []

    def flush():
        if cur is not None and chain:
            out[cur] = np.array([verts[i - 1] for i in chain])

    with open(path, "r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("v "):
                verts.append([float(x) for x in line.split()[1:4]])
            elif line.startswith("o "):
                flush()
                cur = line[2:].strip()
                chain = []
            elif line.startswith("l "):
                idx = [int(x) for x in line.split()[1:]]
                if not chain:
                    chain.extend(idx)
                else:
                    chain.extend(idx[1:] if idx[0] == chain[-1] else idx)
    flush()
    return out
