"""Selection, coloring, slicing and other everyday modifiers."""
from __future__ import annotations

import inspect

import numpy as np

from ..data import COLOR, ID, RADIUS, SELECTION, TRANSPARENCY, TYPE, Frame, scalar_columns
from ..expressions import evaluate
from ..pipeline import Context, Modifier, Param
from ..styles import particle_colors

COLORMAPS = ("viridis", "plasma", "inferno", "magma", "cividis", "turbo", "coolwarm", "RdBu_r",
             "seismic", "jet", "rainbow", "hot", "Blues", "Reds", "Greys", "Spectral")


def _column(frame: Frame, name: str) -> np.ndarray:
    cols = scalar_columns(frame)
    if name in cols:
        return np.asarray(cols[name], dtype=float)
    lower = {k.lower(): k for k in cols}
    if name.lower() in lower:
        return np.asarray(cols[lower[name.lower()]], dtype=float)
    if name in frame.properties and np.asarray(frame.properties[name]).ndim == 2:
        return np.linalg.norm(frame.properties[name], axis=1)
    raise ValueError(f"No property '{name}'. Available: {', '.join(cols)}")


def colormap_rgb(values, cmap, start, end):
    import matplotlib
    cm = matplotlib.colormaps[cmap]
    span = end - start if end != start else 1.0
    t = np.clip((np.asarray(values, dtype=float) - start) / span, 0, 1)
    t = np.nan_to_num(t, nan=0.0)
    return cm(t)[:, :3]


class ColorCoding(Modifier):
    title = "Color coding"
    category = "Coloring"
    description = "Color particles by any property using a colormap, with a legend."
    params = [
        Param("property", "property", "z"),
        Param("colormap", "choice", "viridis", choices=COLORMAPS),
        Param("auto_range", "bool", True, help="Use the min/max of the current frame"),
        Param("robust_range", "bool", False, label="Ignore outliers (1–99%)"),
        Param("start", "float", 0.0, label="Range start"),
        Param("end", "float", 1.0, label="Range end"),
        Param("reverse", "bool", False),
        Param("only_selected", "bool", False),
    ]

    def label(self):
        return f"Color coding: {self.property}"

    def apply(self, frame, ctx):
        v = _column(frame, self.property)
        mask = frame.selection if self.only_selected else np.ones(frame.n, bool)
        vals = v[mask]
        if self.auto_range and len(vals):
            lo, hi = (np.nanpercentile(vals, [1, 99]) if self.robust_range else (np.nanmin(vals), np.nanmax(vals)))
            self.values["start"], self.values["end"] = float(lo), float(hi)
        start, end = self.start, self.end
        cmap = self.colormap + ("_r" if self.reverse and not self.colormap.endswith("_r") else "")
        colors = particle_colors(frame)
        colors[mask] = colormap_rgb(vals, cmap, start, end)
        frame.properties[COLOR] = colors
        if self.only_selected:
            frame.properties[SELECTION] = np.zeros(frame.n, dtype=np.int8)
        frame.attributes["_colorbar"] = {"title": self.property, "cmap": cmap, "start": start, "end": end}
        frame.attributes.pop("_legend", None)
        self.status = f"range {start:.4g} … {end:.4g}"
        return frame


class AssignColor(Modifier):
    title = "Assign color"
    category = "Coloring"
    description = "Give the selected particles a fixed color."
    params = [Param("color", "color", (0.3, 0.5, 1.0)), Param("keep_selection", "bool", False)]

    def apply(self, frame, ctx):
        sel = frame.selection
        colors = particle_colors(frame)
        colors[sel] = np.asarray(self.color, dtype=float)[:3]
        frame.properties[COLOR] = colors
        if not self.keep_selection:
            frame.properties[SELECTION] = np.zeros(frame.n, dtype=np.int8)
        self.status = f"{sel.sum()} particles colored"
        return frame


class ExpressionSelect(Modifier):
    title = "Expression selection"
    category = "Selection"
    description = "Select particles with a Boolean expression, e.g. type == 2 && z > 10."
    params = [Param("expression", "expr", "z > CellOrigin_Z + CellSize_Z / 2")]

    def label(self):
        return f"Select: {self.expression}"

    def apply(self, frame, ctx):
        sel = evaluate(self.expression, frame, ctx.frame_index).astype(bool)
        frame.properties[SELECTION] = sel.astype(np.int8)
        frame.attributes["SelectExpression.num_selected"] = int(sel.sum())
        self.status = f"{sel.sum()} of {frame.n} selected"
        return frame


class SelectType(Modifier):
    title = "Select type"
    category = "Selection"
    description = "Select particles by type number or element name (comma separated)."
    params = [Param("types", "str", "1", help="e.g. 1,3 or Cu,Ni")]

    def apply(self, frame, ctx):
        wanted = [w.strip() for w in str(self.types).split(",") if w.strip()]
        names = {v: k for k, v in frame.type_names.items()}
        ids = set()
        for w in wanted:
            if w.lstrip("-").isdigit():
                ids.add(int(w))
            elif w in names:
                ids.add(names[w])
        types = frame.properties.get(TYPE, np.ones(frame.n, int))
        sel = np.isin(types, list(ids))
        frame.properties[SELECTION] = sel.astype(np.int8)
        self.status = f"{sel.sum()} selected"
        return frame


class InvertSelection(Modifier):
    title = "Invert selection"
    category = "Selection"

    def apply(self, frame, ctx):
        frame.properties[SELECTION] = (~frame.selection).astype(np.int8)
        return frame


class ClearSelection(Modifier):
    title = "Clear selection"
    category = "Selection"

    def apply(self, frame, ctx):
        frame.properties.pop(SELECTION, None)
        return frame


class DeleteSelected(Modifier):
    title = "Delete selected"
    category = "Modification"

    def apply(self, frame, ctx):
        sel = frame.selection
        self.status = f"{sel.sum()} deleted"
        out = frame.subset(~sel)
        out.properties.pop(SELECTION, None)
        return out


class Slice(Modifier):
    title = "Slice"
    category = "Modification"
    description = "Cut away particles on one side of a plane, or keep a slab of given width."
    params = [
        Param("normal", "vector", (1.0, 0.0, 0.0)),
        Param("distance", "float", 0.0, help="Plane position along the normal, from the cell center"),
        Param("width", "float", 0.0, min=0.0, help="0 = half-space; >0 = keep a slab this thick"),
        Param("invert", "bool", False),
        Param("select_only", "bool", False, label="Select instead of delete"),
    ]

    def apply(self, frame, ctx):
        n = np.asarray(self.normal, dtype=float)
        if not np.linalg.norm(n):
            raise ValueError("Normal must be non-zero")
        n = n / np.linalg.norm(n)
        center = frame.origin + frame.cell.sum(axis=0) / 2
        s = (frame.positions - center) @ n - self.distance
        cut = (np.abs(s) > self.width / 2) if self.width > 0 else (s > 0)
        if self.invert:
            cut = ~cut
        if self.select_only:
            frame.properties[SELECTION] = cut.astype(np.int8)
            self.status = f"{cut.sum()} selected"
            return frame
        self.status = f"{cut.sum()} removed"
        return frame.subset(~cut)


class ComputeProperty(Modifier):
    title = "Compute property"
    category = "Modification"
    description = "Create or overwrite a property from an expression, e.g. sqrt(vx^2+vy^2+vz^2)."
    params = [
        Param("output", "str", "my_property"),
        Param("expression", "expr", "x * 0"),
        Param("only_selected", "bool", False),
    ]

    def label(self):
        return f"Compute: {self.output}"

    def apply(self, frame, ctx):
        val = evaluate(self.expression, frame, ctx.frame_index)
        name = self.output.strip()
        target = {"Color.R": (COLOR, 0), "Color.G": (COLOR, 1), "Color.B": (COLOR, 2),
                  "Position.X": ("pos", 0), "Position.Y": ("pos", 1), "Position.Z": ("pos", 2)}.get(name)
        mask = frame.selection if self.only_selected else slice(None)
        if target:
            prop, comp = target
            arr = frame.positions if prop == "pos" else particle_colors(frame)
            arr[mask, comp] = val[mask]
            if prop == COLOR:
                frame.properties[COLOR] = arr
        else:
            if self.only_selected and name in frame.properties:
                arr = np.asarray(frame.properties[name], dtype=float).copy()
                arr[mask] = val[mask]
                val = arr
            frame.properties[name] = val
        self.status = f"min {np.nanmin(val):.4g}, max {np.nanmax(val):.4g}" if val.dtype.kind in "biuf" and len(val) else ""
        return frame


class WrapPeriodic(Modifier):
    title = "Wrap at periodic boundaries"
    category = "Modification"

    def apply(self, frame, ctx):
        f = frame.fractional()
        for d in range(3):
            if frame.pbc[d]:
                f[:, d] %= 1.0
        frame.positions = frame.cartesian(f)
        return frame


class SetCell(Modifier):
    title = "Set cell"
    category = "Modification"
    description = "Fixed orthogonal box for files without one (e.g. XYZ). Length 0 keeps that axis."
    params = [Param("lx", "float", 0.0, min=0), Param("ly", "float", 0.0, min=0), Param("lz", "float", 0.0, min=0),
              Param("centered", "bool", True, help="box from -L/2 to L/2; off: from 0 to L"),
              Param("pbc", "str", "", help="periodic axes, e.g. y or xyz; empty keeps the current setting")]

    def apply(self, frame, ctx):
        o, c = frame.origin.copy(), frame.cell.copy()
        for d, L in enumerate((self.lx, self.ly, self.lz)):
            if L > 0:
                c[d] = 0; c[d, d] = L; o[d] = -L / 2 if self.centered else 0.0
        frame.cell, frame.origin = c, o
        if self.pbc:
            frame.pbc = tuple(a in self.pbc.lower() for a in "xyz")
        return frame


class Replicate(Modifier):
    title = "Replicate"
    category = "Modification"
    description = "Tile the periodic cell n times along each cell vector."
    params = [Param("nx", "int", 2, min=1, max=50), Param("ny", "int", 2, min=1, max=50),
              Param("nz", "int", 1, min=1, max=50), Param("adjust_box", "bool", True)]

    def apply(self, frame, ctx):
        reps = [(i, j, k) for i in range(self.nx) for j in range(self.ny) for k in range(self.nz)]
        shifts = np.array(reps) @ frame.cell
        n = frame.n
        pos = (frame.positions[None, :, :] + shifts[:, None, :]).reshape(-1, 3)
        props = {}
        for k, v in frame.properties.items():
            props[k] = np.concatenate([v] * len(reps))
        if ID in props:
            offset = int(frame.properties[ID].max()) if n else 0
            props[ID] = props[ID] + np.repeat(np.arange(len(reps)), n) * offset
        out = Frame(pos, frame.cell.copy(), frame.origin.copy(), frame.pbc, props, dict(frame.attributes),
                    tables=dict(frame.tables), type_names=dict(frame.type_names))
        if frame.bonds is not None and len(frame.bonds):
            out.bonds = np.concatenate([frame.bonds + r * n for r in range(len(reps))])
        if self.adjust_box:
            out.cell = frame.cell * np.array([[self.nx], [self.ny], [self.nz]])
        return out


class SetRadius(Modifier):
    title = "Set radius"
    category = "Coloring"
    description = "Set the display radius of particles (all or selected)."
    params = [Param("radius", "float", 1.0, min=0.01, step=0.05), Param("only_selected", "bool", False)]

    def apply(self, frame, ctx):
        from ..styles import particle_radii
        r = particle_radii(frame)
        mask = frame.selection if self.only_selected else np.ones(frame.n, bool)
        r[mask] = self.radius
        frame.properties[RADIUS] = r
        return frame


class SetTransparency(Modifier):
    title = "Set transparency"
    category = "Coloring"
    description = "Make particles see-through (0 = opaque, 1 = invisible)."
    params = [Param("transparency", "float", 0.7, min=0.0, max=1.0, step=0.05), Param("only_selected", "bool", True)]

    def apply(self, frame, ctx):
        t = np.asarray(frame.properties.get(TRANSPARENCY, np.zeros(frame.n)), dtype=float).copy()
        mask = frame.selection if self.only_selected else np.ones(frame.n, bool)
        t[mask] = self.transparency
        frame.properties[TRANSPARENCY] = t
        if self.only_selected:
            frame.properties[SELECTION] = np.zeros(frame.n, dtype=np.int8)
        return frame


class Displacement(Modifier):
    title = "Displacement vectors"
    category = "Analysis"
    description = "Displacement of every particle relative to a reference frame (matched by id)."
    params = [Param("reference_frame", "int", 0, min=0), Param("minimum_image", "bool", True)]

    def apply(self, frame, ctx: Context):
        ref = ctx.upstream(int(self.reference_frame))
        cur_pos, ref_pos = frame.positions, ref.positions
        if ID in frame.properties and ID in ref.properties:
            lookup = {int(i): k for k, i in enumerate(ref.properties[ID])}
            idx = np.array([lookup.get(int(i), -1) for i in frame.properties[ID]])
            if (idx < 0).any():
                raise ValueError("Some particles are missing from the reference frame")
            ref_pos = ref_pos[idx]
        elif ref.n != frame.n:
            raise ValueError("Particle count differs from the reference frame and there are no ids")
        d = cur_pos - ref_pos
        if self.minimum_image:
            d = frame.minimum_image(d)
        frame.properties["displacement"] = d
        frame.properties["displacement_magnitude"] = np.linalg.norm(d, axis=1)
        msd = float(np.mean(np.sum(d * d, axis=1))) if frame.n else 0.0
        frame.attributes["MSD"] = msd
        self.status = f"MSD = {msd:.4g}"
        return frame


class Histogram(Modifier):
    title = "Histogram"
    category = "Analysis"
    description = "Histogram of a property, shown in the Plots tab."
    params = [Param("property", "property", "z"), Param("bins", "int", 50, min=2, max=10000),
              Param("only_selected", "bool", False)]

    def apply(self, frame, ctx):
        v = _column(frame, self.property)
        if self.only_selected:
            v = v[frame.selection]
        counts, edges = np.histogram(v[np.isfinite(v)], bins=int(self.bins))
        frame.tables[f"Histogram: {self.property}"] = {
            "x": (edges[:-1] + edges[1:]) / 2, "y": counts, "xlabel": self.property, "ylabel": "Count", "kind": "bar"}
        return frame


class Function(Modifier):
    """Wraps a Python callable so scripts can do pipeline.add(my_function)."""

    title = "Python function"
    category = "Scripting"

    def __init__(self, func=None, **kw):
        super().__init__(**kw)
        self.func = func

    def label(self):
        return f"Python: {getattr(self.func, '__name__', 'function')}"

    def apply(self, frame, ctx):
        nargs = len(inspect.signature(self.func).parameters)
        out = self.func(ctx.frame_index, frame) if nargs >= 2 else self.func(frame)
        return out if isinstance(out, Frame) else frame


class Vectors(Modifier):
    title = "Vectors (arrows)"
    category = "Visual"
    description = ("Arrows from a vector property (velocity, force, displacement, N x 2 or N x 3) "
                   "or from an angle in radians (2D orientation, e.g. theta).")
    params = [Param("property", "property", "velocity"), Param("scale", "float", 1.0, min=0.0, step=0.1),
              Param("width", "float", 0.1, min=0.0, step=0.02, help="shaft radius as a fraction of the arrow length"),
              Param("color_by_particle", "bool", True), Param("color", "color", (0.15, 0.15, 0.15)),
              Param("centered", "bool", False, help="center arrows on particles instead of starting there"),
              Param("unit_length", "bool", False, help="draw every arrow with length = scale")]

    def label(self):
        return f"Vectors: {self.property}"

    def apply(self, frame, ctx):
        v = vector_field(frame, self.property)
        frame.attributes.setdefault("_vectors", [])
        frame.attributes["_vectors"] = frame.attributes["_vectors"] + [dict(self.values)]
        mag = np.linalg.norm(v, axis=1)
        self.status = f"|v| mean {mag.mean():.4g}, max {mag.max():.4g}" if frame.n else ""
        return frame


def vector_field(frame, name):
    # N x 3 vectors from an N x 3 or N x 2 property, or from an angle (radians) in the xy plane.
    if name not in frame.properties:
        cand = {k.lower(): k for k in frame.properties}
        if name.lower() not in cand: raise ValueError(f"No property '{name}'. Available: {', '.join(frame.properties)}")
        name = cand[name.lower()]
    a = np.asarray(frame.properties[name], dtype=float)
    if a.ndim == 1: return np.column_stack([np.cos(a), np.sin(a), np.zeros(len(a))])
    if a.shape[1] == 2: return np.column_stack([a, np.zeros(len(a))])
    return a[:, :3]
