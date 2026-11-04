"""Per-frame analyses for colloids, active matter and other soft systems."""
from __future__ import annotations

import numpy as np

from ..data import COLOR, ID, Frame
from ..neighbors import canonical_system, is2d, pairs_within
from ..pipeline import Modifier, Param
from .basic import vector_field


class StructureFactor(Modifier):
    title = "Structure factor S(k)"
    category = "Analysis"
    description = "Static structure factor from the particle positions (direct method on the reciprocal lattice of the cell)."
    params = [Param("k_max", "float", 10.0, min=0.01, step=0.5), Param("k_min", "float", 0.0, min=0.0, step=0.1),
              Param("bins", "int", 100, min=5, max=5000), Param("only_selected", "bool", False)]

    def apply(self, frame, ctx):
        import freud
        if not all(frame.pbc[:2 if is2d(frame) else 3]): raise ValueError("S(k) needs a periodic cell")
        f = frame.subset(frame.selection) if self.only_selected else frame
        box, pts = canonical_system(f)
        s = freud.diffraction.StaticStructureFactorDirect(bins=int(self.bins), k_max=float(self.k_max), k_min=float(self.k_min))
        s.compute((box, pts), reset=True)
        k, S = np.asarray(s.bin_centers), np.asarray(s.S_k, dtype=float)
        S[k < s.min_valid_k] = np.nan  # below the smallest wave vector the cell supports
        frame.tables["S(k)"] = {"x": k, "y": S, "xlabel": "k", "ylabel": "S(k)", "kind": "line"}
        i = np.nanargmax(S) if np.isfinite(S).any() else 0
        frame.attributes["StructureFactor.peak_k"], frame.attributes["StructureFactor.peak_S"] = float(k[i]), float(S[i])
        self.status = f"first peak k = {k[i]:.4g}, S = {S[i]:.4g}"
        return frame


def orientation_units(frame: Frame, name: str, dim: int) -> np.ndarray:
    u = vector_field(frame, name)[:, :dim]
    nrm = np.linalg.norm(u, axis=1)
    return u / np.where(nrm > 0, nrm, 1)[:, None]


def polar_nematic(u, dim):
    # Global polar order |<u>| and nematic order (largest eigenvalue of Q = (d <uu> - I)/(d - 1)).
    P = float(np.linalg.norm(u.mean(0)))
    Q = (dim * np.einsum("ni,nj->ij", u, u) / len(u) - np.eye(dim)) / (dim - 1)
    return P, float(np.linalg.eigvalsh(Q)[-1])


class OrientationalOrder(Modifier):
    title = "Polar and nematic order"
    category = "Analysis"
    description = ("Global polar order |<u>| and nematic order S of an orientation (vector property or 2D angle), "
                   "and optionally the local polar order within a cutoff.")
    params = [Param("property", "property", "theta"), Param("local_cutoff", "float", 0.0, min=0.0, step=0.5,
                                                           help="0 = global only")]

    def apply(self, frame, ctx):
        dim = 2 if is2d(frame) else 3
        u = orientation_units(frame, self.property, dim)
        P, S = polar_nematic(u, dim)
        frame.attributes["Order.polar"], frame.attributes["Order.nematic"] = P, S
        if self.local_cutoff > 0:
            i, j, _, _ = pairs_within(frame, float(self.local_cutoff))
            acc = u.copy()
            for a in range(dim):
                acc[:, a] += np.bincount(i, u[j, a], frame.n)
            frame.properties["local_polar_order"] = np.linalg.norm(acc, axis=1) / (np.bincount(i, minlength=frame.n) + 1)
        self.status = f"polar {P:.4f}, nematic {S:.4f}"
        return frame


class BondOrderCorrelation(Modifier):
    title = "Bond-orientational correlation g6(r) (2D)"
    category = "Analysis"
    description = "g_k(r) = <Re psi_k(i) psi_k*(j)> over pairs at distance r; its decay distinguishes liquid, hexatic and solid."
    params = [Param("k", "int", 6, min=2, max=12), Param("r_max", "float", 10.0, min=0.1, step=0.5),
              Param("bins", "int", 100, min=5, max=5000), Param("num_neighbors", "int", 6, min=1, max=50)]

    def apply(self, frame, ctx):
        import freud
        if not is2d(frame): raise ValueError("g6(r) needs 2D data")
        box, pts = canonical_system(frame)
        psi = np.asarray(freud.order.Hexatic(k=int(self.k)).compute(
            (box, pts), neighbors={"num_neighbors": int(self.num_neighbors), "exclude_ii": True}).particle_order)
        i, j, d, _ = pairs_within(frame, float(self.r_max))
        edges = np.linspace(0, float(self.r_max), int(self.bins) + 1)
        b = np.clip(np.digitize(d, edges) - 1, 0, int(self.bins) - 1)
        num, cnt = np.bincount(b, (psi[i] * psi[j].conj()).real, int(self.bins)), np.bincount(b, minlength=int(self.bins))
        g = np.where(cnt > 0, num / np.maximum(cnt, 1), np.nan)
        frame.tables[f"g{int(self.k)}(r)"] = {"x": (edges[:-1] + edges[1:]) / 2, "y": g, "xlabel": "r",
                                             "ylabel": f"g{int(self.k)}(r)", "kind": "line"}
        self.status = f"{len(d) // 2} pairs"
        return frame


class Trails(Modifier):
    title = "Particle trails"
    category = "Visual"
    description = "Draws each particle's path over the previous frames (unwrapped across periodic boundaries)."
    params = [Param("length", "int", 10, min=1, max=1000, help="number of previous frames"),
              Param("every", "int", 1, min=1, max=1000), Param("only_selected", "bool", False),
              Param("width", "float", 1.5, min=0.5, max=20, step=0.5)]

    def apply(self, frame, ctx):
        keep = frame.selection if self.only_selected else np.ones(frame.n, bool)
        frames = [f for f in range(ctx.frame_index - int(self.length) * int(self.every), ctx.frame_index, int(self.every)) if f >= 0]
        ids = frame.properties.get(ID)
        path = np.full((len(frames) + 1, frame.n, 3), np.nan)
        path[-1] = frame.positions
        for k in range(len(frames) - 1, -1, -1):  # walk backwards, unwrapping each step against the later frame
            g = ctx.upstream(frames[k])
            if ids is not None and ID in g.properties:
                look = {int(v): n for n, v in enumerate(g.properties[ID])}
                idx = np.array([look.get(int(v), -1) for v in ids])
            else:
                idx = np.arange(frame.n) if g.n == frame.n else -np.ones(frame.n, int)
            ok = idx >= 0
            step = np.full((frame.n, 3), np.nan)
            step[ok] = g.minimum_image(g.positions[idx[ok]] - path[k + 1][ok])
            path[k] = path[k + 1] + step
        frame.attributes["_trails"] = {"path": path[:, keep], "rows": np.nonzero(keep)[0], "width": float(self.width)}
        self.status = f"{len(frames)} previous frames"
        return frame
