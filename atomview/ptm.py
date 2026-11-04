"""Polyhedral template matching (Larsen, Schmidt & Schiotz 2016) through the bundled MIT-licensed PTM library."""
from __future__ import annotations

import ctypes
import os

import numpy as np

from .data import Frame
from .neighbors import nearest_chunks

FLAGS = {"fcc": 1, "hcp": 2, "bcc": 4, "ico": 8, "sc": 16, "dcub": 32, "dhex": 64, "graphene": 128}
NAMES = {0: "Other", 1: "FCC", 2: "HCP", 3: "BCC", 4: "ICO", 5: "SC", 6: "Cubic diamond", 7: "Hex diamond", 8: "Graphene"}
ALLOYS = {0: "None", 1: "Pure", 2: "L10", 3: "L12 (Cu)", 4: "L12 (Au)", 5: "B2", 6: "SiC", 7: "BN"}
MAX_NBRS = 18
_lib = None


def library():
    global _lib
    if _lib is None:
        from . import _ptmlib
        _lib = ctypes.CDLL(_ptmlib.__file__)
        d, i64, i32 = (np.ctypeslib.ndpointer(t, flags="C") for t in (np.float64, np.int64, np.int32))
        _lib.av_ptm.argtypes = [ctypes.c_int64, ctypes.c_int, i64, d, ctypes.c_void_p, ctypes.c_int32, ctypes.c_int,
                                ctypes.c_int, i32, i32, d, d, d, d, d, d, i32]
        _lib.av_ptm.restype = ctypes.c_int
    return _lib


def neighbor_table(frame: Frame, k=MAX_NBRS, cancel=None):
    # Sorted k nearest neighbors and their minimum-image vectors in the frame's own (unrotated) coordinates.
    idx = np.empty((frame.n, k), dtype=np.int64)
    vec = np.zeros((frame.n, k, 3))
    for s, ix, _, _ in nearest_chunks(frame, k, 20000):
        if cancel and cancel():
            from .pipeline import Cancelled
            raise Cancelled()
        m = len(ix)
        idx[s:s + m] = ix
        ok = ix >= 0
        r, c = np.nonzero(ok)
        vec[s + r, c] = frame.minimum_image(frame.positions[ix[ok]] - frame.positions[s + r])
    return idx, vec


def ptm(frame: Frame, structures=("fcc", "hcp", "bcc", "ico"), numbers=None, threads=None, conventional=True, cancel=None):
    """Structure type, RMSD, orientation quaternion (w, x, y, z), interatomic distance, lattice constant,
    elastic deformation gradient F (row-major 3x3) and alloy type for every particle."""
    n = frame.n
    flags = sum(FLAGS[s] for s in structures)
    idx, vec = neighbor_table(frame, MAX_NBRS, cancel)
    out = {"type": np.zeros(n, np.int32), "alloy": np.zeros(n, np.int32), "rmsd": np.zeros(n), "scale": np.zeros(n),
           "q": np.zeros((n, 4)), "F": np.zeros((n, 9)), "interatomic_distance": np.zeros(n),
           "lattice_constant": np.zeros(n), "template": np.zeros(n, np.int32)}
    nums = None if numbers is None else np.ascontiguousarray(numbers, dtype=np.int32)
    library().av_ptm(n, MAX_NBRS, idx, vec, None if nums is None else nums.ctypes.data, flags, int(conventional),
                     threads or os.cpu_count() or 1, out["type"], out["alloy"], out["rmsd"], out["scale"], out["q"], out["F"],
                     out["interatomic_distance"], out["lattice_constant"], out["template"])
    other = out["type"] == 0
    for key in ("rmsd", "scale", "interatomic_distance", "lattice_constant"):
        out[key][other] = np.nan
    out["q"][other] = np.nan
    out["F"][other] = np.nan
    return out
