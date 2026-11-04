"""Neighbor finding on top of freud, for any cell shape and boundary mix.

freud wants an upper-triangular box centred on the origin. We map the frame's
cell to that canonical form (distances are preserved) and pad non-periodic
directions so no particle sits outside the box.
"""
from __future__ import annotations

import freud
import numpy as np

from .data import Frame


def is2d(frame: Frame) -> bool:
    return frame.attributes.get("Dimensions") == 2


def canonical_system(frame: Frame, r_max: float = 0.0):
    """Return (freud.Box, points) for neighbor queries on this frame."""
    cell = frame.cell.copy()
    frac = frame.fractional()
    two = is2d(frame)
    for d in range(2 if two else 3):
        if not frame.pbc[d]:
            length = np.linalg.norm(cell[d])
            fmin, fmax = frac[:, d].min(), frac[:, d].max()
            new_len = (fmax - fmin) * length + 2 * r_max + 2.0
            s = new_len / length
            frac[:, d] = ((frac[:, d] - fmin) * length + r_max + 1.0) / new_len
            cell[d] *= s
        else:
            frac[:, d] %= 1.0
    a, b, c = cell
    lx = np.linalg.norm(a)
    ahat = a / lx
    xy = b @ ahat
    ly = np.sqrt(max(b @ b - xy * xy, 1e-30))
    if two:
        box = freud.box.Box(Lx=lx, Ly=ly, xy=xy / ly, is2D=True)
        box.periodic = [frame.pbc[0], frame.pbc[1], True]  # z is unused in a 2D box but freud requires the flag
        f2 = frac[:, :2] - 0.5
        pts = np.column_stack([f2[:, 0] * lx + f2[:, 1] * xy, f2[:, 1] * ly, np.zeros(len(f2))])
        return box, pts.astype(np.float32)
    xz = c @ ahat
    yz = (b @ c - xy * xz) / ly
    lz = np.sqrt(max(c @ c - xz * xz - yz * yz, 1e-30))
    canon = np.array([[lx, 0, 0], [xy, ly, 0], [xz, yz, lz]])
    box = freud.box.Box.from_matrix(canon.T)
    box.periodic = list(frame.pbc)
    points = ((frac - 0.5) @ canon).astype(np.float32)
    return box, points


def _check_cutoff(box: freud.box.Box, frame: Frame, r: float):
    L = np.array([box.Lx, box.Ly, box.Lz])
    for d in range(3):
        if frame.pbc[d] and r >= L[d] / 2:
            raise ValueError(f"Cutoff {r:g} is too large for the periodic cell (must be < half the cell, "
                             f"{L[d] / 2:.3g}). Use the Replicate modifier first.")


def pairs_within(frame: Frame, cutoff: float):
    """All pairs (i, j) with i != j and distance < cutoff. Each pair appears twice."""
    if frame.n == 0:
        e = np.zeros(0, dtype=np.int64)
        return e, e, np.zeros(0), np.zeros((0, 3))
    box, pts = canonical_system(frame, cutoff)
    _check_cutoff(box, frame, cutoff)
    nl = freud.locality.AABBQuery(box, pts).query(pts, dict(r_max=cutoff, exclude_ii=True)).toNeighborList()
    i = np.asarray(nl.query_point_indices, dtype=np.int64)
    j = np.asarray(nl.point_indices, dtype=np.int64)
    vec = box.wrap(pts[j] - pts[i])
    return i, j, np.asarray(nl.distances, dtype=float), np.asarray(vec, dtype=float)


def nearest_chunks(frame: Frame, k: int, chunk: int = 20_000):
    """k nearest neighbors, sorted by distance, yielded per block of particles so memory stays O(chunk).

    Yields (start, idx (m, k), vec (m, k, 3), dist (m, k)) for particles start..start+m; rows with fewer
    than k neighbors are padded with -1 / nan.
    """
    n = frame.n
    box, pts = canonical_system(frame)
    q = freud.locality.AABBQuery(box, pts)
    kk = min(k + 1, n)  # +1: the particle itself, removed by global index below
    for s in range(0, n, chunk):
        e = min(n, s + chunk)
        nl = q.query(pts[s:e], dict(num_neighbors=kk)).toNeighborList()
        qi = np.asarray(nl.query_point_indices, dtype=np.int64)
        pj = np.asarray(nl.point_indices, dtype=np.int64)
        dist = np.asarray(nl.distances, dtype=float)
        keep = pj != qi + s
        qi, pj, dist = qi[keep], pj[keep], dist[keep]
        order = np.lexsort((dist, qi))
        qi, pj, dist = qi[order], pj[order], dist[order]
        m = e - s
        starts = np.concatenate([[0], np.cumsum(np.bincount(qi, minlength=m))[:-1]])
        rank = np.arange(len(qi)) - starts[qi]
        ok = rank < k
        idx = -np.ones((m, k), dtype=np.int64)
        d = np.full((m, k), np.nan)
        idx[qi[ok], rank[ok]] = pj[ok]
        d[qi[ok], rank[ok]] = dist[ok]
        vec = np.full((m, k, 3), np.nan)
        rows = np.nonzero(idx >= 0)
        vec[rows] = box.wrap(pts[idx[rows]] - pts[s + rows[0]])
        yield s, idx, vec, d


def nearest(frame: Frame, k: int):
    """k nearest neighbors of every particle: (indices (N, k), vectors (N, k, 3), distances (N, k))."""
    parts = list(nearest_chunks(frame, k, max(frame.n, 1)))
    if not parts:
        return np.zeros((0, k), dtype=np.int64), np.zeros((0, k, 3)), np.zeros((0, k))
    return parts[0][1], parts[0][2], parts[0][3]
