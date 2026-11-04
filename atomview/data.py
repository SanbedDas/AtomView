"""The particle dataset that flows through a pipeline.

A Frame holds one snapshot: positions, the simulation cell, per-particle
properties (any number of named arrays), global attributes, optional bonds and
data tables (for plots such as an RDF or a histogram).
"""
from __future__ import annotations

import copy
import re
from dataclasses import dataclass, field

import numpy as np

# Standard property names used across readers, modifiers and the renderer.
ID = "id"
TYPE = "type"
COLOR = "color"
RADIUS = "radius"
SELECTION = "selection"
TRANSPARENCY = "transparency"


def sanitize(name: str) -> str:
    """Turn a property name like 'f_1[2]' or 'Position.X' into an identifier."""
    s = re.sub(r"[^0-9a-zA-Z_]", "_", name).strip("_")
    if not s:
        s = "prop"
    if s[0].isdigit():
        s = "p_" + s
    return s


@dataclass
class Frame:
    positions: np.ndarray
    cell: np.ndarray = field(default_factory=lambda: np.eye(3))  # rows are a, b, c
    origin: np.ndarray = field(default_factory=lambda: np.zeros(3))
    pbc: tuple = (True, True, True)
    properties: dict = field(default_factory=dict)
    attributes: dict = field(default_factory=dict)
    bonds: np.ndarray | None = None  # (M, 2) particle indices
    tables: dict = field(default_factory=dict)  # name -> {"x":…, "y":…, "xlabel":…, "ylabel":…}
    type_names: dict = field(default_factory=dict)  # int type -> name (element)

    def __post_init__(self):
        self.positions = np.ascontiguousarray(self.positions, dtype=float).reshape(-1, 3)
        self.cell = np.asarray(self.cell, dtype=float).reshape(3, 3)
        self.origin = np.asarray(self.origin, dtype=float).reshape(3)
        self.pbc = tuple(bool(p) for p in self.pbc)

    # ---- convenience -------------------------------------------------------
    @property
    def n(self) -> int:
        return len(self.positions)

    def __len__(self):
        return self.n

    def __getitem__(self, name):
        if name in ("position", "positions", "Position"):
            return self.positions
        return self.properties[name]

    def __setitem__(self, name, value):
        if name in ("position", "positions", "Position"):
            self.positions = np.asarray(value, dtype=float).reshape(-1, 3)
            return
        arr = np.asarray(value)
        if arr.ndim == 0:
            arr = np.full(self.n, arr)
        if len(arr) != self.n:
            raise ValueError(f"Property '{name}' has {len(arr)} values but there are {self.n} particles")
        self.properties[name] = arr

    def __contains__(self, name):
        return name in self.properties

    def get(self, name, default=None):
        return self.properties.get(name, default)

    def copy(self) -> "Frame":
        return copy.deepcopy(self)

    def subset(self, keep: np.ndarray) -> "Frame":
        """Return a new frame with only the particles where `keep` is true."""
        keep = np.asarray(keep, dtype=bool)
        f = Frame(
            positions=self.positions[keep],
            cell=self.cell.copy(),
            origin=self.origin.copy(),
            pbc=self.pbc,
            properties={k: v[keep] for k, v in self.properties.items()},
            attributes=dict(self.attributes),
            tables=dict(self.tables),
            type_names=dict(self.type_names),
        )
        if self.bonds is not None and len(self.bonds):
            new_index = -np.ones(self.n, dtype=np.int64)
            new_index[keep] = np.arange(keep.sum())
            b = new_index[self.bonds]
            f.bonds = b[(b >= 0).all(axis=1)]
        return f

    @property
    def selection(self) -> np.ndarray:
        s = self.properties.get(SELECTION)
        if s is None:
            return np.zeros(self.n, dtype=bool)
        return np.asarray(s) != 0

    def fractional(self, positions=None) -> np.ndarray:
        p = self.positions if positions is None else positions
        return (p - self.origin) @ np.linalg.inv(self.cell)

    def cartesian(self, frac) -> np.ndarray:
        return frac @ self.cell + self.origin

    def minimum_image(self, vectors) -> np.ndarray:
        """Apply the minimum image convention to displacement vectors."""
        inv = np.linalg.inv(self.cell)
        f = vectors @ inv
        for d in range(3):
            if self.pbc[d]:
                f[:, d] -= np.round(f[:, d])
        return f @ self.cell

    def volume(self) -> float:
        return abs(np.linalg.det(self.cell))

    def summary(self) -> str:
        lengths = np.linalg.norm(self.cell, axis=1)
        props = ", ".join(self.properties)
        return (f"{self.n} particles, cell {lengths[0]:.3g} × {lengths[1]:.3g} × {lengths[2]:.3g}, "
                f"properties: {props}")


def scalar_columns(frame: Frame) -> dict:
    """Every per-particle quantity as a flat 1D array, keyed by an identifier.

    Vector properties are split into components (color_R, displacement_X, …).
    Positions appear as x, y, z. Used for expressions, color coding and tables.
    """
    cols = {"x": frame.positions[:, 0], "y": frame.positions[:, 1], "z": frame.positions[:, 2]}
    for name, arr in frame.properties.items():
        key = sanitize(name)
        arr = np.asarray(arr)
        if arr.ndim == 1:
            cols[key] = arr
        elif arr.ndim == 2:
            comps = "RGBA" if name == COLOR else "XYZW"
            for i in range(arr.shape[1]):
                suffix = comps[i] if arr.shape[1] <= 4 else str(i)
                cols[f"{key}_{suffix}"] = arr[:, i]
    return cols
