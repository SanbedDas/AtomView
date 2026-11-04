"""Default per-type colors and radii, and structure-type colors."""
from __future__ import annotations

import numpy as np

from .data import COLOR, RADIUS, TYPE, Frame

# Pleasant, distinguishable defaults for numeric LAMMPS types (1-based).
TYPE_PALETTE = np.array([
    [0.89, 0.35, 0.29], [0.30, 0.56, 0.86], [0.98, 0.78, 0.25], [0.40, 0.73, 0.42],
    [0.62, 0.43, 0.80], [0.95, 0.55, 0.20], [0.35, 0.80, 0.80], [0.85, 0.45, 0.65],
    [0.60, 0.60, 0.60], [0.55, 0.40, 0.30],
])

STRUCTURE_NAMES = {0: "Other", 1: "FCC", 2: "HCP", 3: "BCC", 4: "ICO", 5: "SC", 6: "Cubic diamond", 7: "Hex diamond", 8: "Graphene"}
STRUCTURE_COLORS = {
    0: (0.95, 0.95, 0.95),
    1: (0.40, 1.00, 0.40),
    2: (1.00, 0.40, 0.40),
    3: (0.40, 0.40, 1.00),
    4: (0.95, 0.75, 0.10),
    5: (0.63, 0.43, 0.90),
    6: (0.07, 0.75, 0.85),
    7: (1.00, 0.60, 0.15),
    8: (0.55, 0.55, 0.55),
}


def _element(name):
    from ase.data import atomic_numbers
    if not name:
        return None
    name = name.strip()
    cand = name[:2].capitalize() if len(name) > 1 and name[1].islower() else name[:1].upper()
    return atomic_numbers.get(name.capitalize()) or atomic_numbers.get(cand)


def type_color(t: int, name: str = "") -> np.ndarray:
    z = _element(name)
    if z:
        from ase.data.colors import jmol_colors
        return np.array(jmol_colors[z])
    return TYPE_PALETTE[(int(t) - 1) % len(TYPE_PALETTE)]


def type_radius(t: int, name: str = "", fallback: float = 0.5) -> float:
    z = _element(name)
    if z:
        from ase.data import covalent_radii
        return float(covalent_radii[z])
    return fallback


def auto_radius(frame: Frame) -> float:
    """A radius that makes spheres nearly touch for a dense solid."""
    if frame.n == 0 or frame.volume() <= 0:
        return 0.5
    if frame.attributes.get("Dimensions") == 2:
        spacing = (np.linalg.norm(np.cross(frame.cell[0], frame.cell[1])) / frame.n) ** 0.5
    else:
        spacing = (frame.volume() / frame.n) ** (1 / 3)
    return float(np.clip(0.5 * spacing, 0.2, 3.0))


def particle_colors(frame: Frame) -> np.ndarray:
    if COLOR in frame.properties:
        c = np.asarray(frame.properties[COLOR], dtype=float)
        return c[:, :3] if c.ndim == 2 else np.tile(c[:, None], (1, 3))
    types = frame.properties.get(TYPE)
    if types is None:
        return np.tile(TYPE_PALETTE[0], (frame.n, 1))
    out = np.empty((frame.n, 3))
    for t in np.unique(types):
        out[types == t] = type_color(t, frame.type_names.get(int(t), ""))
    return out


def particle_radii(frame: Frame) -> np.ndarray:
    if RADIUS in frame.properties:
        return np.asarray(frame.properties[RADIUS], dtype=float)
    types = frame.properties.get(TYPE)
    fallback = auto_radius(frame)
    if types is None:
        return np.full(frame.n, fallback)
    out = np.empty(frame.n)
    for t in np.unique(types):
        out[types == t] = type_radius(t, frame.type_names.get(int(t), ""), fallback)
    return out


def type_legend(frame: Frame):
    types = frame.properties.get(TYPE)
    if types is None or COLOR in frame.properties:
        return []
    return [(frame.type_names.get(int(t), f"Type {int(t)}"), type_color(t, frame.type_names.get(int(t), "")))
            for t in np.unique(types)]
