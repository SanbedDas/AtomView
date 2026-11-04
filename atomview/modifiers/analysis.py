"""Structural analysis modifiers: coordination/RDF, adaptive CNA, centrosymmetry,
clusters, bonds, Steinhardt order parameters and Voronoi analysis."""
from __future__ import annotations

import numpy as np

from ..data import COLOR, SELECTION, Frame
from ..neighbors import canonical_system, is2d, nearest, nearest_chunks, pairs_within
from ..pipeline import Cancelled, Modifier, Param
from ..styles import STRUCTURE_COLORS, STRUCTURE_NAMES


class Coordination(Modifier):
    title = "Coordination analysis"
    category = "Analysis"
    description = "Neighbor count within a cutoff, and the radial distribution function g(r)."
    params = [Param("cutoff", "float", 3.2, min=0.01, step=0.1), Param("bins", "int", 200, min=10, max=5000)]

    def apply(self, frame, ctx):
        import freud
        i, _, _, _ = pairs_within(frame, self.cutoff)
        frame.properties["coordination"] = np.bincount(i, minlength=frame.n)
        box, pts = canonical_system(frame, self.cutoff)
        if all(frame.pbc[:2 if is2d(frame) else 3]):
            rdf = freud.density.RDF(bins=int(self.bins), r_max=float(self.cutoff))
            rdf.compute(system=(box, pts), reset=True)
            frame.tables["Radial distribution function"] = {
                "x": np.asarray(rdf.bin_centers), "y": np.asarray(rdf.rdf), "xlabel": "r", "ylabel": "g(r)", "kind": "line"}
        mean = frame.properties["coordination"].mean() if frame.n else 0
        self.status = f"mean coordination {mean:.3f}"
        return frame


def _adjacency(vec, rc):
    diff = vec[:, :, None, :] - vec[:, None, :, :]
    A = np.einsum("mabx,mabx->mab", diff, diff) <= (rc ** 2)[:, None, None]
    k = vec.shape[1]
    A[:, np.arange(k), np.arange(k)] = False
    return A.astype(np.int32)


def _signatures(A):
    ncn = A.sum(2)
    AA = A @ A  # AA[a, b] = number of common neighbors of a and b
    deg = A * AA  # degree of b inside the common-neighbor graph of bond a
    nb = deg.sum(2) // 2
    return ncn, nb, deg.max(2)


def _two_triangles(A, rows):
    """True where the 6-ring common-neighbor graph of bond (m, a) is actually two triangles."""
    S = A[rows]  # (M, k, k)
    M, k, _ = S.shape
    out = np.zeros((M, k), dtype=bool)
    for m in range(M):
        Am = S[m]
        for a in range(k):
            c = np.nonzero(Am[a])[0]
            if len(c) != 6:
                continue
            sub = Am[np.ix_(c, c)]
            out[m, a] = np.trace(sub @ sub @ sub) > 0
    return out


def adaptive_cna(frame: Frame, cancel=None) -> np.ndarray:
    """Adaptive common neighbor analysis (Stukowski 2012). Returns structure ids
    0 Other, 1 FCC, 2 HCP, 3 BCC, 4 ICO."""
    n = frame.n
    result = np.zeros(n, dtype=np.int32)
    if n < 15:
        return result
    for s, idx, vec, dist in nearest_chunks(frame, 14, 20000):
        if cancel and cancel(): raise Cancelled()
        res = result[s:s + len(idx)]  # view: writes go to result
        rows = np.arange(len(idx))
        valid12, valid14 = (idx[:, :12] >= 0).all(axis=1), (idx >= 0).all(axis=1)
        # FCC / HCP / ICO from the 12 nearest neighbors
        ok = rows[valid12]
        if len(ok):
            rc = 1.2071068 * np.mean(dist[ok, :12], axis=1)
            A = _adjacency(vec[ok, :12], rc)
            ncn, nb, maxdeg = _signatures(A)
            s421 = (ncn == 4) & (nb == 2) & (maxdeg == 1)
            s422 = (ncn == 4) & (nb == 2) & (maxdeg == 2)
            s555 = (ncn == 5) & (nb == 5) & (maxdeg == 2)
            n421, n422, n555 = s421.sum(1), s422.sum(1), s555.sum(1)
            res[ok[n421 == 12]] = 1
            res[ok[(n421 == 6) & (n422 == 6)]] = 2
            res[ok[n555 == 12]] = 4
        # BCC from the 14 nearest neighbors, for the rest
        left = rows[valid14 & (res == 0)]
        if len(left):
            scale = (np.sum(dist[left, :8], axis=1) / np.sqrt(0.75) + np.sum(dist[left, 8:14], axis=1)) / 14
            rc = 1.2071068 * scale
            A = _adjacency(vec[left], rc)
            ncn, nb, maxdeg = _signatures(A)
            s666 = (ncn == 6) & (nb == 6) & (maxdeg == 2)
            s444 = (ncn == 4) & (nb == 4) & (maxdeg == 2)
            cand = (s666.sum(1) == 8) & (s444.sum(1) == 6)
            if cand.any():
                tri = _two_triangles(A, np.nonzero(cand)[0])
                bad = (tri & s666[cand]).any(axis=1)
                res[left[np.nonzero(cand)[0][~bad]]] = 3
    return result


class CommonNeighborAnalysis(Modifier):
    title = "Common neighbor analysis"
    category = "Analysis"
    description = "Adaptive CNA: identifies FCC, HCP, BCC and icosahedral atoms (no cutoff needed)."
    params = [Param("color_by_type", "bool", True)]

    def apply(self, frame, ctx):
        if is2d(frame): raise ValueError("Common neighbor analysis is for 3D crystals; use Hexatic order for 2D")
        st = adaptive_cna(frame, ctx.cancelled if ctx else None)
        frame.properties["structure_type"] = st
        counts = np.bincount(st, minlength=5)
        for k in range(5):  # CNA types are the first five structure ids
            frame.attributes[f"CNA.counts.{STRUCTURE_NAMES[k]}"] = int(counts[k])
        if self.color_by_type:
            frame.properties[COLOR] = np.array([STRUCTURE_COLORS[int(t)] for t in range(5)])[st]
            frame.attributes["_legend"] = [(f"{STRUCTURE_NAMES[k]} ({counts[k]})", STRUCTURE_COLORS[k])
                                           for k in range(5) if counts[k]]
            frame.attributes.pop("_colorbar", None)
        n = max(frame.n, 1)
        self.status = ", ".join(f"{STRUCTURE_NAMES[k]} {100 * counts[k] / n:.1f}%" for k in range(5) if counts[k])
        return frame


class CentroSymmetry(Modifier):
    title = "Centrosymmetry parameter"
    category = "Analysis"
    description = "Kelchner centrosymmetry: ~0 in perfect FCC/BCC, large at defects and surfaces."
    params = [Param("num_neighbors", "int", 12, min=2, max=40, step=2)]

    def apply(self, frame, ctx):
        k = int(self.num_neighbors) // 2 * 2
        a, b = np.triu_indices(k, 1)
        csp = np.empty(frame.n)
        for s, idx, vec, _ in nearest_chunks(frame, k, 20000):
            sums = vec[:, a, :] + vec[:, b, :]
            d2 = np.einsum("mpx,mpx->mp", sums, sums)
            d2 = np.where(np.isnan(d2), np.inf, d2)
            d2.sort(axis=1)
            csp[s:s + len(idx)] = d2[:, : k // 2].sum(axis=1)
        csp[~np.isfinite(csp)] = np.nan
        frame.properties["centrosymmetry"] = csp
        self.status = f"median {np.nanmedian(csp):.3g}" if frame.n else ""
        return frame


class ClusterAnalysis(Modifier):
    title = "Cluster analysis"
    category = "Analysis"
    description = "Groups particles into clusters of neighbors within a cutoff."
    params = [Param("cutoff", "float", 3.2, min=0.01, step=0.1), Param("only_selected", "bool", False),
              Param("sort_by_size", "bool", True), Param("color_clusters", "bool", False)]

    def apply(self, frame, ctx):
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import connected_components
        mask = frame.selection if self.only_selected else np.ones(frame.n, bool)
        sub = frame.subset(mask)
        i, j, _, _ = pairs_within(sub, self.cutoff)
        g = coo_matrix((np.ones(len(i)), (i, j)), shape=(sub.n, sub.n))
        ncl, labels = connected_components(g, directed=False)
        if self.sort_by_size and ncl:
            sizes = np.bincount(labels)
            rank = np.empty(ncl, dtype=np.int64)
            rank[np.argsort(-sizes, kind="stable")] = np.arange(ncl)
            labels = rank[labels]
        cluster = np.zeros(frame.n, dtype=np.int64)
        cluster[mask] = labels + 1
        frame.properties["cluster"] = cluster
        sizes = np.bincount(labels) if ncl else np.zeros(0)
        frame.attributes["ClusterAnalysis.cluster_count"] = int(ncl)
        frame.attributes["ClusterAnalysis.largest_size"] = int(sizes.max()) if ncl else 0
        frame.tables["Cluster sizes"] = {"x": np.arange(1, ncl + 1), "y": np.sort(sizes)[::-1] if ncl else sizes,
                                         "xlabel": "Cluster (by size)", "ylabel": "Particles", "kind": "bar"}
        if self.color_clusters and ncl:
            rng = np.random.default_rng(1)
            palette = rng.uniform(0.2, 0.95, size=(ncl + 1, 3))
            palette[0] = (0.8, 0.8, 0.8)
            frame.properties[COLOR] = palette[cluster]
        if self.only_selected:
            frame.properties[SELECTION] = np.zeros(frame.n, dtype=np.int8)
        self.status = f"{ncl} clusters, largest {frame.attributes['ClusterAnalysis.largest_size']}"
        return frame


class CreateBonds(Modifier):
    title = "Create bonds"
    category = "Visual"
    description = "Bonds between particles closer than a cutoff (optionally per type pair, e.g. 1-2:2.1)."
    params = [Param("cutoff", "float", 3.0, min=0.01, step=0.1),
              Param("pair_cutoffs", "str", "", help="Optional: '1-2:2.1, O-H:1.2' overrides the cutoff"),
              Param("bond_radius", "float", 0.15, min=0.0, step=0.05)]

    def apply(self, frame, ctx):
        from ..data import TYPE
        pairs = {}
        names = {v: k for k, v in frame.type_names.items()}
        for item in str(self.pair_cutoffs).split(","):
            if ":" not in item:
                continue
            p, c = item.split(":")
            a, b = [w.strip() for w in p.split("-")]
            ta = int(a) if a.lstrip("-").isdigit() else names.get(a)
            tb = int(b) if b.lstrip("-").isdigit() else names.get(b)
            if ta is not None and tb is not None:
                pairs[(ta, tb)] = pairs[(tb, ta)] = float(c)
        rmax = max([self.cutoff] + list(pairs.values()))
        i, j, d, _ = pairs_within(frame, rmax)
        keep = i < j
        i, j, d = i[keep], j[keep], d[keep]
        if pairs:
            types = frame.properties.get(TYPE, np.ones(frame.n, int))
            ti, tj = types[i], types[j]
            cut = np.array([pairs.get((int(a), int(b)), 0.0) for a, b in zip(ti, tj)])
            ok = d < cut
        else:
            ok = d < self.cutoff
        frame.bonds = np.column_stack([i[ok], j[ok]])
        frame.attributes["_bond_radius"] = float(self.bond_radius)
        self.status = f"{len(frame.bonds)} bonds"
        return frame


class Steinhardt(Modifier):
    title = "Steinhardt order parameter"
    category = "Analysis"
    description = "Bond-orientational order q_l (e.g. q6 to tell solid from liquid)."
    params = [Param("l", "int", 6, min=2, max=12), Param("num_neighbors", "int", 12, min=1, max=100),
              Param("average", "bool", False, label="Averaged (Lechner–Dellago)")]

    def apply(self, frame, ctx):
        import freud
        if is2d(frame): raise ValueError("Steinhardt q_l is defined in 3D; use Hexatic order for 2D")
        box, pts = canonical_system(frame)
        st = freud.order.Steinhardt(int(self.l), average=bool(self.average))
        st.compute((box, pts), neighbors={"num_neighbors": int(self.num_neighbors), "exclude_ii": True})
        name = f"{'avg_' if self.average else ''}q{int(self.l)}"
        frame.properties[name] = np.asarray(st.particle_order, dtype=float)
        self.status = f"mean {name} = {np.nanmean(frame.properties[name]):.3f}"
        return frame


class VoronoiAnalysis(Modifier):
    title = "Voronoi analysis"
    category = "Analysis"
    description = "Voronoi cell volume and number of faces per particle (periodic cells)."
    params = [Param("min_face_area", "float", 0.0, min=0.0, step=0.01, help="Ignore tiny faces when counting")]

    def apply(self, frame, ctx):
        import freud
        if not all(frame.pbc[:2 if is2d(frame) else 3]):
            raise ValueError("Voronoi analysis needs a periodic cell")
        box, pts = canonical_system(frame)
        vor = freud.locality.Voronoi()
        vor.compute((box, pts))
        frame.properties["voronoi_volume"] = np.asarray(vor.volumes, dtype=float)
        nl = vor.nlist
        w = np.asarray(nl.weights)
        qi = np.asarray(nl.query_point_indices)
        frame.properties["voronoi_coordination"] = np.bincount(qi[w > self.min_face_area], minlength=frame.n)
        self.status = f"total volume {frame.properties['voronoi_volume'].sum():.4g}"
        return frame


class Hexatic(Modifier):
    title = "Hexatic order (2D)"
    category = "Analysis"
    description = "Bond-orientational order psi_k of 2D systems: |psi6| ~ 1 in a hexagonal crystal, ~0 in a liquid."
    params = [Param("k", "int", 6, min=2, max=12), Param("num_neighbors", "int", 6, min=1, max=50)]

    def apply(self, frame, ctx):
        import freud
        if not is2d(frame): raise ValueError("Hexatic order needs 2D data (all particles in one z plane)")
        box, pts = canonical_system(frame)
        h = freud.order.Hexatic(k=int(self.k)).compute((box, pts), neighbors={"num_neighbors": int(self.num_neighbors), "exclude_ii": True})
        psi = np.asarray(h.particle_order)
        k = int(self.k)
        frame.properties[f"psi{k}_abs"], frame.properties[f"psi{k}_phase"] = np.abs(psi), np.angle(psi)
        frame.attributes[f"Hexatic.global_psi{k}"] = float(np.abs(psi.mean()))
        frame.attributes[f"Hexatic.mean_abs_psi{k}"] = float(np.abs(psi).mean())
        self.status = f"<|psi{k}|> = {np.abs(psi).mean():.3f}, |<psi{k}>| = {np.abs(psi.mean()):.3f}"
        return frame


def atomic_strain(ref: Frame, cur: Frame, cutoff, dim=3):
    # Per-atom deformation gradient F from neighbors within cutoff in ref (Shimizu, Ogata & Li 2007), with
    # Green-Lagrangian strain E = (F^T F - I)/2 and the non-affine residual D2min (Falk & Langer 1998).
    # Rows of ref and cur are the same particles. Returns (F (N,3,3), E (N,3,3), D2min (N,), valid (N,)).
    n = ref.n
    i, j, _, _ = pairs_within(ref, cutoff)
    d0 = ref.minimum_image(ref.positions[j] - ref.positions[i])[:, :dim]
    d = cur.minimum_image(cur.positions[j] - cur.positions[i])[:, :dim]
    V, W = np.zeros((n, dim, dim)), np.zeros((n, dim, dim))
    for a in range(dim):
        for b in range(dim):
            V[:, a, b] = np.bincount(i, d0[:, a] * d0[:, b], n)
            W[:, a, b] = np.bincount(i, d0[:, a] * d[:, b], n)
    valid = np.abs(np.linalg.det(V)) > 1e-8 * np.maximum(np.einsum("nii->n", V), 1e-300) ** dim
    J = np.full((n, dim, dim), np.nan)  # d = d0 @ J for each pair, so F = J^T
    J[valid] = np.linalg.solve(V[valid], W[valid])
    r = d - np.einsum("ma,mab->mb", d0, J[i])
    d2 = np.bincount(i, (r ** 2).sum(1), n)
    d2[~valid] = np.nan
    F = np.tile(np.eye(3), (n, 1, 1))
    F[:, :dim, :dim] = np.transpose(J, (0, 2, 1))
    F[~valid] = np.nan
    E = 0.5 * (np.einsum("nki,nkj->nij", F, F) - np.eye(3))
    return F, E, d2, valid


def shear_strain(E, dim=3):
    # von Mises shear strain invariant (2D uses only in-plane components).
    if dim == 2:
        return np.sqrt(E[:, 0, 1] ** 2 + (E[:, 0, 0] - E[:, 1, 1]) ** 2 / 2)
    return np.sqrt(E[:, 0, 1] ** 2 + E[:, 0, 2] ** 2 + E[:, 1, 2] ** 2 +
                   ((E[:, 1, 1] - E[:, 2, 2]) ** 2 + (E[:, 0, 0] - E[:, 2, 2]) ** 2 + (E[:, 0, 0] - E[:, 1, 1]) ** 2) / 6)


class AtomicStrain(Modifier):
    title = "Atomic strain"
    category = "Analysis"
    description = ("Local deformation relative to a reference frame: von Mises shear strain, volumetric strain, "
                   "Green-Lagrangian strain tensor and non-affine D2min. Particles are matched by id.")
    params = [Param("cutoff", "float", 3.0, min=0.01, step=0.1, help="neighbor cutoff in the reference frame"),
              Param("reference_frame", "int", 0, min=0), Param("output_tensor", "bool", False, label="Output strain tensor"),
              Param("output_gradient", "bool", False, label="Output deformation gradient")]

    def apply(self, frame, ctx):
        from ..data import ID, SELECTION
        ref = ctx.upstream(int(self.reference_frame))
        if ID in frame.properties and ID in ref.properties:
            look = {int(v): k for k, v in enumerate(ref.properties[ID])}
            idx = np.array([look.get(int(v), -1) for v in frame.properties[ID]])
            if (idx < 0).any(): raise ValueError(f"{int((idx < 0).sum())} particles are missing from the reference frame")
            ref = _reorder(ref, idx)
        elif ref.n != frame.n:
            raise ValueError("Particle count differs from the reference frame and there are no ids")
        dim = 2 if is2d(frame) else 3
        F, E, d2, valid = atomic_strain(ref, frame, float(self.cutoff), dim)
        frame.properties["shear_strain"] = shear_strain(E, dim)
        frame.properties["volumetric_strain"] = np.trace(E[:, :dim, :dim], axis1=1, axis2=2) / dim
        frame.properties["D2min"] = d2
        if self.output_tensor:
            frame.properties["strain_tensor"] = E[:, [0, 1, 2, 0, 0, 1], [0, 1, 2, 1, 2, 2]]  # xx yy zz xy xz yz
        if self.output_gradient:
            frame.properties["deformation_gradient"] = F.reshape(-1, 9)
        frame.properties[SELECTION] = (~valid).astype(np.int8)
        frame.attributes["AtomicStrain.invalid_count"] = int((~valid).sum())
        self.status = (f"mean shear strain {np.nanmean(frame.properties['shear_strain']):.4g}; "
                       f"{int((~valid).sum())} particles with too few neighbors selected")
        return frame


def _reorder(f: Frame, idx):
    g = f.subset(np.ones(f.n, bool))
    g.positions = f.positions[idx]
    g.properties = {k: np.asarray(v)[idx] for k, v in f.properties.items()}
    return g


class PolyhedralTemplateMatching(Modifier):
    title = "Polyhedral template matching"
    category = "Analysis"
    description = ("PTM (Larsen et al. 2016): structure type with RMSD, local lattice orientation, interatomic distance, "
                   "elastic deformation gradient and chemical ordering. More robust than CNA at high temperature and strain.")
    params = [Param("fcc", "bool", True, label="FCC"), Param("hcp", "bool", True, label="HCP"),
              Param("bcc", "bool", True, label="BCC"), Param("ico", "bool", True, label="Icosahedral"),
              Param("sc", "bool", False, label="Simple cubic"), Param("dcub", "bool", False, label="Cubic diamond"),
              Param("dhex", "bool", False, label="Hexagonal diamond"), Param("graphene", "bool", False, label="Graphene"),
              Param("rmsd_cutoff", "float", 0.1, min=0.0, step=0.01, help="matches with larger RMSD count as Other; 0 = keep all"),
              Param("output_orientation", "bool", True), Param("output_deformation_gradient", "bool", False),
              Param("output_alloy", "bool", False, label="Output ordering (uses particle types)"),
              Param("color_by_type", "bool", True)]

    def apply(self, frame, ctx):
        from ..ptm import ALLOYS, ptm
        from ..data import TYPE
        sel = [s for s in ("fcc", "hcp", "bcc", "ico", "sc", "dcub", "dhex", "graphene") if self.values[s]]
        if not sel: raise ValueError("Select at least one structure")
        nums = frame.properties.get(TYPE) if self.output_alloy else None
        r = ptm(frame, sel, nums, cancel=ctx.cancelled if ctx else None)
        st, rmsd = r["type"].astype(np.int32), r["rmsd"]
        if self.rmsd_cutoff > 0:
            st[rmsd > self.rmsd_cutoff] = 0
        frame.properties["structure_type"] = st
        frame.properties["ptm_rmsd"] = rmsd
        frame.properties["interatomic_distance"] = r["interatomic_distance"]
        if self.output_orientation:
            frame.properties["orientation"] = r["q"]  # quaternion (w, x, y, z)
        if self.output_deformation_gradient:
            frame.properties["elastic_deformation_gradient"] = r["F"]
        if self.output_alloy:
            frame.properties["ordering_type"] = r["alloy"]
            frame.attributes.update({f"PTM.ordering.{v}": int((r["alloy"] == k).sum()) for k, v in ALLOYS.items() if (r["alloy"] == k).any()})
        counts = np.bincount(st, minlength=9)
        for k, name in STRUCTURE_NAMES.items():
            frame.attributes[f"PTM.counts.{name.replace(' ', '_')}"] = int(counts[k])
        ok = np.isfinite(rmsd)
        if ok.any():
            h, e = np.histogram(rmsd[ok], bins=60)
            frame.tables["PTM RMSD"] = {"x": (e[:-1] + e[1:]) / 2, "y": h, "xlabel": "RMSD", "ylabel": "Count", "kind": "bar"}
        if self.color_by_type:
            frame.properties[COLOR] = np.array([STRUCTURE_COLORS[k] for k in range(9)])[st]
            frame.attributes["_legend"] = [(f"{STRUCTURE_NAMES[k]} ({counts[k]})", STRUCTURE_COLORS[k]) for k in range(9) if counts[k]]
            frame.attributes.pop("_colorbar", None)
        n = max(frame.n, 1)
        self.status = ", ".join(f"{STRUCTURE_NAMES[k]} {100 * counts[k] / n:.1f}%" for k in range(9) if counts[k])
        return frame
