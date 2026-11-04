"""Reading and writing particle files.

Native readers: LAMMPS dump (text, .gz, triclinic, scaled/unwrapped coords),
LAMMPS data, XYZ / extended XYZ, VASP POSCAR/CONTCAR. Everything else ASE can
read (CIF, PDB, XDATCAR, OUTCAR, GRO, XSF, CFG, …) goes through ASE.

Trajectories are indexed once (byte offsets of each frame) and frames are
parsed on demand, so multi-gigabyte dumps open instantly.
"""
from __future__ import annotations

import glob
import gzip
import os
import re
from collections import OrderedDict

import numpy as np

from .data import ID, TYPE, Frame


def _open(path, mode="rb"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def _natural_key(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s.lower())]


# --------------------------------------------------------------------------
# Trajectory source
# --------------------------------------------------------------------------
class Trajectory:
    """A sequence of frames read lazily from one or more files."""

    def __init__(self, path: str, topology: str | None = None):
        self.path = str(path)
        self.topology = topology
        if any(c in self.path for c in "*?["):
            files = sorted(glob.glob(self.path), key=_natural_key)
            if not files:
                raise FileNotFoundError(f"No files match {self.path}")
        else:
            if not os.path.exists(self.path):
                raise FileNotFoundError(self.path)
            files = [self.path]
        self.files = files
        self.format = detect_format(files[0])
        self._index = []  # list of (file, offset or frame number)
        for f in files:
            for loc in (_index_file(f, self.format) if self.format != "chemfiles" else range(chem_reader(f, topology).n)):
                self._index.append((f, loc))
        if not self._index:
            raise ValueError(f"No frames found in {self.path}")
        self._cache: OrderedDict[int, Frame] = OrderedDict()
        self._ase_cache = {}

    def __len__(self):
        return len(self._index)

    @property
    def n_frames(self):
        return len(self._index)

    def __getitem__(self, i: int) -> Frame:
        return self.get(i)

    def get(self, i: int) -> Frame:
        i = int(i) % len(self._index)
        if i in self._cache:
            self._cache.move_to_end(i)
            return self._cache[i]
        file, loc = self._index[i]
        frame = (chem_reader(file, self.topology).frame(loc) if self.format == "chemfiles"
                 else _read_frame(file, self.format, loc, self._ase_cache))
        frame.attributes.setdefault("SourceFrame", i)
        frame.attributes.setdefault("SourceFile", os.path.basename(file))
        mark_2d(frame)
        self._cache[i] = frame
        while len(self._cache) > 8:
            self._cache.popitem(last=False)
        return frame


def load(path: str, topology: str | None = None) -> Trajectory:
    return Trajectory(path, topology)


def mark_2d(f: Frame):
    # A frame with every particle in one z plane is treated as a 2D system (z non-periodic).
    if f.attributes.get("Dimensions") is None and f.n > 2 and np.ptp(f.positions[:, 2]) == 0:
        f.attributes["Dimensions"] = 2
    if f.attributes.get("Dimensions") == 2:
        f.pbc = (f.pbc[0], f.pbc[1], False)
    return f


def detect_format(path: str) -> str:
    name = os.path.basename(path).lower()
    base = name[:-3] if name.endswith(".gz") else name
    if base.endswith(".xyz") or base.endswith(".extxyz"):
        return "xyz"
    if base.endswith(".mat"):
        return "mat"
    if base.endswith(CHEMFILES_EXT):
        return "chemfiles"
    if base.startswith(("poscar", "contcar")) or base.endswith((".poscar", ".vasp")):
        return "poscar"
    if base.endswith((".data", ".lmp", ".lammps-data")) or base.startswith("data."):
        if _peek(path).startswith("ITEM:"):
            return "lammps_dump"
        return "lammps_data"
    head = _peek(path)
    if head.startswith("ITEM:"):
        return "lammps_dump"
    if any(base.endswith(e) for e in (".dump", ".lammpstrj", ".lammps", ".atom")):
        return "lammps_dump"
    lines = head.splitlines()
    if lines and lines[0].strip().isdigit():
        return "xyz"
    return "ase"


def _peek(path, n=4096) -> str:
    try:
        with _open(path) as fh:
            return fh.read(n).decode("utf-8", "replace").lstrip()
    except OSError:
        return ""


def _index_file(path, fmt):
    if fmt == "lammps_dump":
        return _index_pattern(path, b"ITEM: TIMESTEP")
    if fmt == "xyz":
        return _index_xyz(path)
    if fmt == "mat":
        return list(range(len(mat_reader(path))))
    if fmt == "ase":
        from ase.io import read
        try:
            n = len(read(path, index=":"))
        except Exception as e:
            raise ValueError(f"Could not read {path}: {e}") from e
        return list(range(n))
    return [0]


def _index_pattern(path, marker: bytes):
    offsets = []
    with _open(path) as fh:
        pos = 0
        for line in fh:
            if line.startswith(marker):
                offsets.append(pos)
            pos += len(line)
    return offsets


def _index_xyz(path):
    offsets = []
    with _open(path) as fh:
        pos = 0
        while True:
            line = fh.readline()
            if not line:
                break
            stripped = line.strip()
            if not stripped:
                pos += len(line)
                continue
            try:
                n = int(stripped.split()[0])
            except ValueError:
                raise ValueError(f"Bad XYZ atom count line in {path}: {stripped[:50]!r}")
            offsets.append(pos)
            pos += len(line)
            for _ in range(n + 1):
                pos += len(fh.readline())
    return offsets


def _read_frame(path, fmt, loc, ase_cache) -> Frame:
    if fmt == "lammps_dump":
        with _open(path) as fh:
            fh.seek(loc)
            return _parse_lammps_dump(fh)
    if fmt == "xyz":
        with _open(path) as fh:
            fh.seek(loc)
            return _parse_xyz(fh)
    if fmt == "lammps_data":
        return read_lammps_data(path)
    if fmt == "poscar":
        return read_poscar(path)
    if fmt == "mat":
        return mat_reader(path).frame(loc)
    if path not in ase_cache:
        from ase.io import read
        ase_cache.clear()
        ase_cache[path] = read(path, index=":")
    return from_ase(ase_cache[path][loc])


def _columns_from_tokens(lines, ncol):
    """Parse whitespace columns; numeric columns become float/int arrays."""
    tokens = b" ".join(lines).split()
    n = len(lines)
    if len(tokens) != n * ncol:
        raise ValueError(f"Expected {ncol} columns per line")
    try:
        arr = np.array(tokens, dtype=float).reshape(n, ncol)
        return [arr[:, i] for i in range(ncol)]
    except ValueError:
        t = np.array(tokens).reshape(n, ncol)
        cols = []
        for i in range(ncol):
            try:
                cols.append(t[:, i].astype(float))
            except ValueError:
                cols.append(t[:, i].astype(str))
        return cols


def _as_int_if_integral(a):
    if a.dtype.kind == "f" and len(a) and np.all(a == np.round(a)) and np.abs(a).max() < 2**62:
        return a.astype(np.int64)
    return a


def _lammps_cell(bounds, tilts):
    (xlo_b, xhi_b), (ylo_b, yhi_b), (zlo_b, zhi_b) = bounds
    xy, xz, yz = tilts
    xlo = xlo_b - min(0.0, xy, xz, xy + xz)
    xhi = xhi_b - max(0.0, xy, xz, xy + xz)
    ylo = ylo_b - min(0.0, yz)
    yhi = yhi_b - max(0.0, yz)
    cell = np.array([[xhi - xlo, 0, 0], [xy, yhi - ylo, 0], [xz, yz, zhi_b - zlo_b]], dtype=float)
    return cell, np.array([xlo, ylo, zlo_b])


def _parse_lammps_dump(fh) -> Frame:
    timestep = None
    n = None
    cell = np.eye(3)
    origin = np.zeros(3)
    pbc = (True, True, True)
    while True:
        line = fh.readline()
        if not line:
            raise ValueError("Unexpected end of LAMMPS dump")
        if line.startswith(b"ITEM: TIMESTEP"):
            timestep = int(fh.readline().split()[0])
        elif line.startswith(b"ITEM: NUMBER OF ATOMS"):
            n = int(fh.readline().split()[0])
        elif line.startswith(b"ITEM: BOX BOUNDS"):
            words = line.decode().split()[3:]
            triclinic = "xy" in words
            flags = [w for w in words if w not in ("xy", "xz", "yz", "abc", "origin")]
            if len(flags) >= 3:
                pbc = tuple(f.startswith("p") for f in flags[:3])
            rows = [list(map(float, fh.readline().split())) for _ in range(3)]
            if "abc" in words:  # general triclinic (newer LAMMPS)
                cell = np.array([r[:3] for r in rows])
                origin = np.array([r[3] for r in rows])
            else:
                tilts = (rows[0][2], rows[1][2], rows[2][2]) if triclinic else (0.0, 0.0, 0.0)
                cell, origin = _lammps_cell([r[:2] for r in rows], tilts)
        elif line.startswith(b"ITEM: ATOMS"):
            names = line.decode().split()[2:]
            lines = [fh.readline() for _ in range(n)]
            cols = _columns_from_tokens(lines, len(names))
            break
    data = dict(zip(names, cols))
    frac_sets = [("xs", "ys", "zs"), ("xsu", "ysu", "zsu")]
    cart_sets = [("x", "y", "z"), ("xu", "yu", "zu")]
    pos = None
    used = set()
    for s in cart_sets:
        if all(c in data for c in s):
            pos = np.column_stack([data[c] for c in s]).astype(float)
            used |= set(s)
            break
    if pos is None:
        for s in frac_sets:
            if all(c in data for c in s):
                frac = np.column_stack([data[c] for c in s]).astype(float)
                pos = frac @ cell + origin
                used |= set(s)
                break
    if pos is None:
        raise ValueError("Dump has no x/y/z, xu/yu/zu or xs/ys/zs columns")
    props = {}
    type_names = {}
    for name, col in data.items():
        if name in used:
            continue
        if name == "element":
            continue
        key = {"id": ID, "type": TYPE, "mol": "molecule", "q": "charge"}.get(name, name)
        props[key] = _as_int_if_integral(col) if key in (ID, TYPE, "molecule") or name.startswith(("ix", "iy", "iz")) else col
    if "element" in data:
        elements = data["element"].astype(str)
        if TYPE not in props:
            uniq = list(dict.fromkeys(elements))
            lookup = {e: i + 1 for i, e in enumerate(uniq)}
            props[TYPE] = np.array([lookup[e] for e in elements])
        for t, e in zip(props[TYPE], elements):
            type_names.setdefault(int(t), e)
    # Combine vx/vy/vz style triples into vectors for convenience.
    for base, comps in (("velocity", ("vx", "vy", "vz")), ("force", ("fx", "fy", "fz"))):
        if all(c in props for c in comps):
            props[base] = np.column_stack([props.pop(c) for c in comps])
    if ID in props:  # LAMMPS dumps are often unsorted; sort by id for stable frames
        order = np.argsort(props[ID], kind="stable")
        if not np.all(order == np.arange(len(order))):
            pos = pos[order]
            props = {k: v[order] for k, v in props.items()}
    return Frame(pos, cell, origin, pbc, props, {"Timestep": timestep}, type_names=type_names)


# --------------------------------------------------------------------------
# XYZ / extended XYZ
# --------------------------------------------------------------------------
_KV = re.compile(r'(\w+)=("[^"]*"|\S+)')


def _parse_xyz(fh) -> Frame:
    n = int(fh.readline().split()[0])
    comment = fh.readline().decode("utf-8", "replace")
    lines = [fh.readline() for _ in range(n)]
    info = {k: v.strip('"') for k, v in _KV.findall(comment)}
    cell = None
    pbc = (False, False, False)
    if "Lattice" in info:
        cell = np.array(info["Lattice"].split(), dtype=float).reshape(3, 3)
        pbc = (True, True, True)
    if "pbc" in info:
        pbc = tuple(w.upper().startswith("T") for w in info["pbc"].split())
    props = {}
    type_names = {}
    if "Properties" in info:
        spec = info["Properties"].split(":")
        fields = [(spec[i], spec[i + 1], int(spec[i + 2])) for i in range(0, len(spec), 3)]
        ncol = sum(f[2] for f in fields)
        cols = _columns_from_tokens(lines, ncol)
        c = 0
        pos = None
        species = None
        for name, kind, count in fields:
            block = cols[c:c + count]
            c += count
            arr = np.column_stack(block) if count > 1 else block[0]
            lname = name.lower()
            if lname in ("pos", "positions"):
                pos = np.asarray(arr, dtype=float)
            elif lname in ("species", "element", "symbols"):
                species = np.asarray(arr).astype(str)
            elif kind == "L":
                props[name] = np.asarray(arr).astype(str) == "T"
            elif kind == "I":
                props[name] = np.asarray(arr).astype(np.int64)
            elif kind == "R":
                key = "velocity" if lname in ("velo", "vel", "velocities") else ("force" if lname == "forces" else name)
                props[key] = np.asarray(arr, dtype=float)
            else:
                props[name] = np.asarray(arr).astype(str)
    else:
        ncol = len(lines[0].split()) if lines else 4
        cols = _columns_from_tokens(lines, ncol)
        species = np.asarray(cols[0]).astype(str)
        pos = np.column_stack(cols[1:4]).astype(float)
        for i, col in enumerate(cols[4:]):
            props[f"col{i + 5}"] = col
    if species is not None:
        if all(s.lstrip("-").isdigit() for s in species[:100]):
            props[TYPE] = species.astype(np.int64)
        else:
            uniq = list(dict.fromkeys(species))
            lookup = {e: i + 1 for i, e in enumerate(uniq)}
            props[TYPE] = np.array([lookup[e] for e in species], dtype=np.int64)
            type_names = {i: e for e, i in lookup.items()}
    origin = np.zeros(3)
    if "Origin" in info:
        origin = np.array(info["Origin"].split(), dtype=float)
    if cell is None:
        lo, hi = pos.min(axis=0), pos.max(axis=0)
        pad = 1.0
        cell = np.diag(np.maximum(hi - lo + 2 * pad, 1.0))
        origin = lo - pad
    attrs = {}
    for k, v in info.items():
        if k in ("Lattice", "Properties", "pbc", "Origin"):
            continue
        try:
            attrs[k] = float(v)
        except ValueError:
            attrs[k] = v
    if "Timestep" not in attrs and "step" in attrs:
        attrs["Timestep"] = attrs["step"]
    return Frame(pos, cell, origin, pbc, props, attrs, type_names=type_names)


# --------------------------------------------------------------------------
# LAMMPS data
# --------------------------------------------------------------------------
_ATOM_STYLES = {
    "atomic": ["id", "type", "x", "y", "z"],
    "charge": ["id", "type", "q", "x", "y", "z"],
    "full": ["id", "mol", "type", "q", "x", "y", "z"],
    "molecular": ["id", "mol", "type", "x", "y", "z"],
    "bond": ["id", "mol", "type", "x", "y", "z"],
    "angle": ["id", "mol", "type", "x", "y", "z"],
    "sphere": ["id", "type", "diameter", "density", "x", "y", "z"],
}


def read_lammps_data(path) -> Frame:
    with _open(path) as fh:
        lines = fh.read().decode("utf-8", "replace").splitlines()
    n = 0
    bounds = [[0, 1], [0, 1], [0, 1]]
    tilts = [0.0, 0.0, 0.0]
    sections = {}
    styles = {}
    i = 1
    header_done = False
    current = None
    while i < len(lines):
        raw = lines[i]
        line = raw.split("#")[0].strip()
        if not header_done:
            if line.endswith(" atoms"):
                n = int(line.split()[0])
            elif line.endswith("xlo xhi"):
                bounds[0] = list(map(float, line.split()[:2]))
            elif line.endswith("ylo yhi"):
                bounds[1] = list(map(float, line.split()[:2]))
            elif line.endswith("zlo zhi"):
                bounds[2] = list(map(float, line.split()[:2]))
            elif line.endswith("xy xz yz"):
                tilts = list(map(float, line.split()[:3]))
        if line and line[0].isalpha() and not line[0].isdigit():
            name = line.split()[0]
            if name in ("Atoms", "Masses", "Velocities", "Bonds"):
                header_done = True
                current = name
                sections[current] = []
                if "#" in raw:
                    styles[current] = raw.split("#")[1].strip()
                i += 1
                continue
            elif header_done:
                current = None
        elif line and current is not None:
            sections[current].append(raw.strip() if current == "Masses" else line)
        i += 1
    rows = sections.get("Atoms", [])[:n] if n else sections.get("Atoms", [])
    ncol = len(rows[0].split())
    style = styles.get("Atoms", "").split()[0] if styles.get("Atoms") else None
    if style not in _ATOM_STYLES:
        style = {5: "atomic", 8: "atomic", 7: "full", 10: "full", 6: "charge", 9: "charge"}.get(ncol, "atomic")
    names = _ATOM_STYLES[style]
    cols = _columns_from_tokens([r.encode() for r in rows], ncol)
    data = dict(zip(names, cols[:len(names)]))
    if ncol >= len(names) + 3:
        for k, c in zip(("ix", "iy", "iz"), cols[len(names):len(names) + 3]):
            data[k] = c
    cell, origin = _lammps_cell(bounds, tilts)
    pos = np.column_stack([data["x"], data["y"], data["z"]]).astype(float)
    props = {ID: data["id"].astype(np.int64), TYPE: data["type"].astype(np.int64)}
    if "mol" in data:
        props["molecule"] = data["mol"].astype(np.int64)
    if "q" in data:
        props["charge"] = data["q"].astype(float)
    for k in ("ix", "iy", "iz"):
        if k in data:
            props[k] = data[k].astype(np.int64)
    # Masses lines often carry the element as a comment: "1 63.546 # Cu"
    masses, type_names = {}, {}
    for r in sections.get("Masses", []):
        body, _, comment = r.partition("#")
        w = body.split()
        masses[int(w[0])] = float(w[1])
        if comment.strip():
            type_names[int(w[0])] = comment.split()[0]
    if masses:
        props["mass"] = np.array([masses.get(int(t), 0.0) for t in props[TYPE]])
    order = np.argsort(props[ID], kind="stable")
    pos = pos[order]
    props = {k: v[order] for k, v in props.items()}
    frame = Frame(pos, cell, origin, (True, True, True), props, {"Timestep": 0}, type_names=type_names)
    if sections.get("Bonds"):
        b = np.array([list(map(int, r.split()[:4])) for r in sections["Bonds"]])
        id_to_index = {int(v): i for i, v in enumerate(props[ID])}
        frame.bonds = np.array([[id_to_index[a], id_to_index[c]] for a, c in b[:, 2:4]], dtype=np.int64)
    return frame


# --------------------------------------------------------------------------
# POSCAR
# --------------------------------------------------------------------------
def read_poscar(path) -> Frame:
    with _open(path) as fh:
        lines = [l for l in fh.read().decode("utf-8", "replace").splitlines()]
    scale = float(lines[1].split()[0])
    lattice = np.array([list(map(float, lines[i].split()[:3])) for i in (2, 3, 4)])
    if scale < 0:
        lattice *= (-scale / abs(np.linalg.det(lattice))) ** (1 / 3)
    else:
        lattice *= scale
    words = lines[5].split()
    if all(w.isdigit() for w in words):
        species = [f"T{i + 1}" for i in range(len(words))]
        counts = list(map(int, words))
        k = 6
    else:
        species = [w.split("/")[0].split("_")[0] for w in words]
        counts = list(map(int, lines[6].split()))
        k = 7
    selective = lines[k].strip().lower().startswith("s")
    if selective:
        k += 1
    direct = lines[k].strip().lower()[0] in "d"
    k += 1
    n = sum(counts)
    coords = np.array([list(map(float, lines[k + i].split()[:3])) for i in range(n)])
    pos = coords @ lattice if direct else coords * (scale if scale > 0 else 1)
    types = np.repeat(np.arange(1, len(counts) + 1), counts)
    return Frame(pos, lattice, np.zeros(3), (True, True, True),
                 {TYPE: types}, {"Timestep": 0, "Comment": lines[0].strip()},
                 type_names={i + 1: s for i, s in enumerate(species)})


# --------------------------------------------------------------------------
# MATLAB .mat files (v4/v5/v7 via SciPy, v7.3/HDF5 via h5py)
# --------------------------------------------------------------------------
_POS_KEYS = ("position", "positions", "pos", "coords", "coordinates", "coord", "xyz", "r", "atoms", "positions_xyz")
_TYPE_KEYS = ("type", "types", "atomtype", "atomtypes", "atom_type", "atom_types", "species", "elements",
              "element", "symbols", "particletype")
_CELL_KEYS = ("cell", "box", "lattice", "h", "boxmatrix", "simbox", "supercell", "lattice_vectors")
_SIZE_KEYS = ("boxsize", "box_size", "boxlength", "box_length", "lengths", "l")
_BOUNDS_KEYS = ("bounds", "boxbounds", "box_bounds", "limits")
_ORIGIN_KEYS = ("origin", "boxlo", "box_lo", "cellorigin")
_TIME_KEYS = ("timestep", "timesteps", "step", "steps")
_ID_KEYS = ("id", "ids", "atomid", "atom_id", "particleidentifier")


class _H5Lazy:
    # HDF5 dataset seen in MATLAB order (reversed axes); slicing reads only the requested hyperslab.
    def __init__(self, ds): self.ds, self.shape, self.ndim, self.dtype, self.size = ds, ds.shape[::-1], ds.ndim, ds.dtype, ds.size

    def __getitem__(self, idx):
        idx = (idx if isinstance(idx, tuple) else (idx,)) + (slice(None),) * self.ndim
        return np.asarray(self.ds[tuple(reversed(idx[:self.ndim]))]).T

    def __array__(self, dtype=None, copy=None): return np.asarray(self.ds[()], dtype=dtype).T


def _load_mat73(path, lazy=False, min_lazy=200_000):
    import h5py

    def convert(f, obj):
        if isinstance(obj, h5py.Group):
            return {k: convert(f, v) for k, v in obj.items() if not k.startswith("#")}
        cls = obj.attrs.get("MATLAB_class", b"")
        cls = cls.decode() if isinstance(cls, bytes) else str(cls)
        if lazy and cls not in ("char", "cell", "logical") and obj.dtype.kind in "iuf" and obj.ndim >= 2 and obj.size >= min_lazy:
            return _H5Lazy(obj)
        data = obj[()]
        if cls == "char":
            return "".join(chr(c) for c in np.asarray(data).ravel())
        if cls == "cell" or (isinstance(data, np.ndarray) and data.dtype == object):
            items = [convert(f, f[ref]) for ref in np.asarray(data).T.ravel()]
            return np.array(items) if all(isinstance(i, str) for i in items) else items
        a = np.asarray(data)
        if cls == "logical":
            a = a.astype(bool)
        return a.T if a.ndim >= 2 else a  # MATLAB stores column-major

    f = h5py.File(path, "r")
    d = {k: convert(f, v) for k, v in f.items() if not k.startswith("#")}
    if not any(isinstance(v, _H5Lazy) for v in _walk(d)):
        f.close()
        f = None
    return d, f


def _walk(d):
    for v in d.values():
        yield from (_walk(v) if isinstance(v, dict) else (v,))


def _load_mat_vars(path, lazy=False):
    # Returns (variables, open h5py file or None). v7.3 files are HDF5; large arrays stay on disk when lazy.
    from scipy.io import loadmat
    try:
        d = loadmat(path, simplify_cells=True)
        d, h5 = {k: v for k, v in d.items() if not k.startswith("__")}, None
    except (NotImplementedError, ValueError):
        d, h5 = _load_mat73(path, lazy)
    # A file holding a single struct (e.g. `data.pos`, `data.type`): use its fields.
    while len(d) == 1 and isinstance(next(iter(d.values())), dict):
        d = next(iter(d.values()))
    return d, h5


def _arr(v):
    return v if isinstance(v, (_H5Lazy, np.ndarray)) else np.asarray(v)


def _positions_frames(pos):
    # Returns (frame count, dimensions, f -> N x 3 positions); 2D input gets z=0. Works on lazy arrays.
    pos = _arr(pos)
    sh = pos.shape
    for d in (3, 2):
        get = None
        if len(sh) == 2 and sh[1] == d: nf, get = 1, lambda f: pos[:, :]
        elif len(sh) == 2 and sh[0] == d and sh[1] > 3: nf, get = 1, lambda f: pos[:, :].T
        elif len(sh) == 3 and sh[1] == d: nf, get = sh[2], lambda f: pos[:, :, f]  # N x d x F
        elif len(sh) == 3 and sh[2] == d: nf, get = sh[0], lambda f: pos[f]  # F x N x d
        elif len(sh) == 3 and sh[0] == d: nf, get = sh[2], lambda f: pos[:, :, f].T  # d x N x F
        if get is not None:
            g = (lambda f: (lambda p: np.column_stack([p, np.zeros(len(p))]))(np.asarray(get(f), dtype=float))) if d == 2 \
                else (lambda f: np.asarray(get(f), dtype=float))
            return nf, d, g
    raise ValueError(f"Positions must be N x 3 or N x 2 (optionally x frames), got {sh}")


def _cell_2d(c, kind):
    # Promotes a 2D cell (2x2 matrix, 2 lengths, or 2x2 lo/hi bounds) to 3D with a unit z axis.
    c = np.asarray(c, dtype=float)
    if kind == "bounds" and c.shape == (2, 2): return np.diag([*(c[:, 1] - c[:, 0]), 1.0]), np.r_[c[:, 0], 0.0]
    if c.shape == (2, 2): return np.block([[c, np.zeros((2, 1))], [np.zeros((1, 2)), np.ones((1, 1))]]), None
    if c.size == 2: return np.diag([*c.ravel(), 1.0]), None
    return None


def _cell_for(value, f, nframes):
    """(cell, origin) from a cell/box variable, possibly per frame."""
    c = np.asarray(value, dtype=float)
    if c.ndim == 3 and c.shape[:2] == (3, 3):
        c = c[:, :, min(f, c.shape[2] - 1)]
    elif c.ndim == 3 and c.shape[1:] == (3, 3):
        c = c[min(f, c.shape[0] - 1)]
    if c.shape == (3, 3):
        return c, None
    if c.size == 3:
        return np.diag(c.ravel()), None
    if c.shape in ((3, 2), (2, 3)):
        b = c if c.shape == (3, 2) else c.T
        return np.diag(b[:, 1] - b[:, 0]), b[:, 0]
    if c.ndim == 2 and 3 in c.shape and nframes in c.shape:  # 3 x F box lengths
        lengths = c[:, f] if c.shape[0] == 3 else c[f]
        return np.diag(lengths), None
    raise ValueError(f"Cannot interpret cell with shape {c.shape}")


class MatReader:
    """Particles from a .mat file, one frame at a time. Variable names are matched case-insensitively.
    For v7.3 (HDF5) files large arrays are read per frame, so memory does not grow with trajectory length."""

    def __init__(self, path, lazy=True):
        d, self.h5 = _load_mat_vars(path, lazy)
        low = {k.lower(): k for k in d}

        def find(keys):
            for k in keys:
                if k in low:
                    return low[k], d[low[k]]
            return None, None

        used = set()
        key, pos = find(_POS_KEYS)
        if pos is None and all(c in low for c in ("x", "y", "z")):
            pos = np.stack([np.asarray(d[low[c]], dtype=float) for c in "xyz"], axis=1)
            used |= {low["x"], low["y"], low["z"]}
        if pos is None:  # last resort: the first numeric N x 3 array
            for k, v in d.items():
                a = _arr(v) if not isinstance(v, (dict, list, str)) else None
                if a is not None and a.dtype.kind in "iuf" and len(a.shape) in (2, 3) and 3 in a.shape[1:] and a.shape[0] > 3:
                    key, pos = k, a
                    break
        if pos is None:
            raise ValueError(f"No particle positions found in {os.path.basename(path)}. Variables: {', '.join(d)}. "
                             "Name the N x 3 array 'Position' (or pos, coords, xyz).")
        used.add(key)
        self.pos_array = _arr(pos)
        self.nf, self.dim, self.pos = _positions_frames(self.pos_array)
        n = self.n = len(self.pos(0))

        tkey, types = find(_TYPE_KEYS)
        self.type_names = {}
        if types is not None:
            used.add(tkey)
            t = np.asarray(types)
            if t.dtype.kind in "USO":
                names = [str(s).strip() for s in t.ravel()]
                lookup = {e: i + 1 for i, e in enumerate(dict.fromkeys(names))}
                types = np.array([lookup[e] for e in names], dtype=np.int64)
                self.type_names = {i: e for e, i in lookup.items()}
            else:
                types = t.reshape(-1).astype(np.int64) if t.size == n else t.astype(np.int64)
        self.types = types
        self.cell_val, self.cell_kind = None, None
        for kind, keys in (("cell", _CELL_KEYS), ("size", _SIZE_KEYS), ("bounds", _BOUNDS_KEYS)):
            ck, cv = find(keys)
            if cv is not None and not isinstance(cv, (dict, str)):
                self.cell_val, self.cell_kind = np.asarray(cv), kind
                used.add(ck)
                break
        for name, keys in (("origin_val", _ORIGIN_KEYS), ("pbc_val", ("pbc", "periodic")), ("steps", _TIME_KEYS)):
            k, v = find(keys)
            setattr(self, name, None if v is None else np.asarray(v))
            if k:
                used.add(k)
        ikey, ids = find(_ID_KEYS)
        self.ids = None
        if ids is not None and np.asarray(ids).size == n:
            self.ids = np.asarray(ids).reshape(-1).astype(np.int64)
            used.add(ikey)
        self.attr_arrays = {k[5:]: np.asarray(v).reshape(-1) for k, v in d.items() if k.startswith("attr_") and not isinstance(v, (dict, list, str))}
        names = d.get("type_names")
        if names is not None and not self.type_names and self.types is not None:
            names = [str(x) for x in np.asarray(names, dtype=object).ravel()]
            self.type_names = {i + 1: nm for i, nm in enumerate(names) if nm}
        self.props = {}
        for k, v in d.items():
            if k in used or k.startswith(("attr_", "table_")) or k == "type_names" or isinstance(v, (dict, list, str)):
                continue
            a = _arr(v)
            if a.dtype.kind not in "biuf" or len(a.shape) == 0 or a.shape[0] != n or n <= 3:
                continue
            self.props[k] = a

    def __len__(self):
        return self.nf

    def frame(self, f) -> Frame:
        p, nf, n, dim = self.pos(f), self.nf, self.n, self.dim
        props = {}
        if self.types is not None:
            props[TYPE] = self.types if self.types.ndim == 1 else self.types[:, min(f, self.types.shape[1] - 1)]
        if self.ids is not None:
            props[ID] = self.ids
        for k, a in self.props.items():
            if len(a.shape) == 3 and a.shape[2] == nf:  # per-frame N x k x F
                props[k] = np.asarray(a[:, :, f])
            elif len(a.shape) == 2 and nf > 1 and a.shape[1] == nf:  # per-frame N x F scalar
                props[k] = np.asarray(a[:, f])
            else:
                props[k] = np.asarray(a)
            if props[k].ndim == 2 and props[k].shape[1] == 1:  # N x 1 column (v7.3 files) -> vector, as v5 files give
                props[k] = props[k][:, 0]
        origin = np.zeros(3)
        c2 = _cell_2d(self.cell_val, self.cell_kind) if (dim == 2 and self.cell_val is not None) else None
        if c2 is not None:
            cell, o = c2
            pbc = (True, True, False)
            if o is not None:
                origin = o
        elif self.cell_val is not None:
            cell, o = _cell_for(self.cell_val, f, nf)
            pbc = (True, True, True)
            if o is not None:
                origin = o
        else:
            lo, hi = p.min(axis=0), p.max(axis=0)
            cell, origin, pbc = np.diag(np.maximum(hi - lo + 2.0, 1.0)), lo - 1.0, (False, False, False)
        if self.origin_val is not None:
            ov = self.origin_val.astype(float).reshape(-1)[:3]
            origin = np.r_[ov, np.zeros(3 - len(ov))]
        if self.pbc_val is not None:
            pv = self.pbc_val.reshape(-1)
            pbc = tuple(bool(x) for x in (pv if pv.size == 3 else (*pv[:2], False) if pv.size == 2 else np.repeat(pv[:1], 3)))
        step = 0
        if self.steps is not None:
            s = self.steps.reshape(-1)
            step = int(s[min(f, len(s) - 1)]) if s.size else 0
        attrs = {"Timestep": step, **({"Dimensions": 2} if dim == 2 else {})}
        for k, arr in self.attr_arrays.items():
            if len(arr) == nf:
                attrs[k] = float(arr[f])
        if dim == 2:
            pbc = (pbc[0], pbc[1], False)
        return Frame(p, cell, origin, pbc, props, attrs, type_names=dict(self.type_names))


_MAT_READERS: dict = {}


def mat_reader(path) -> MatReader:
    key = (os.path.abspath(path), os.path.getmtime(path))
    if key not in _MAT_READERS:
        _MAT_READERS.clear()
        _MAT_READERS[key] = MatReader(path)
    return _MAT_READERS[key]


def read_mat(path) -> list[Frame]:
    """All frames of a .mat file as a list (eager)."""
    r = MatReader(path, lazy=False)
    return [r.frame(f) for f in range(len(r))]


# --------------------------------------------------------------------------
# MD formats through chemfiles (BSD): XTC, TRR, DCD, NetCDF, TNG, GRO, PDB, MOL2, TPR, ...
# --------------------------------------------------------------------------
CHEMFILES_EXT = (".xtc", ".trr", ".dcd", ".nc", ".ncrst", ".tng", ".gro", ".pdb", ".ent", ".mol2", ".tpr", ".arc",
                 ".trj", ".mmtf", ".psf", ".sdf", ".cml", ".cssr", ".molden", ".mmcif", ".smi")
_NO_TOPOLOGY = (".xtc", ".trr", ".dcd", ".nc", ".tng", ".trj")
_TOPOLOGY_EXT = (".gro", ".pdb", ".psf", ".tpr", ".mol2", ".data")


def sibling_topology(path):
    # A same-named .gro/.pdb/.psf/.tpr/.mol2 next to a trajectory that carries no atom names.
    if not path.lower().endswith(_NO_TOPOLOGY): return None
    stem = os.path.splitext(path)[0]
    return next((stem + e for e in _TOPOLOGY_EXT if os.path.exists(stem + e)), None)


class ChemReader:
    def __init__(self, path, topology=None):
        import chemfiles as cf
        self.t = cf.Trajectory(path)
        self.topology = topology or sibling_topology(path)
        if self.topology:
            self.t.set_topology(self.topology)
        self.n = self.t.nsteps
        self.cf = cf

    def frame(self, i) -> Frame:
        f = self.t.read_step(i)
        pos = np.array(f.positions, dtype=float)
        n = len(pos)
        names = [a.name or a.type for a in f.atoms]
        lookup = {e: k + 1 for k, e in enumerate(dict.fromkeys(names))}
        props = {ID: np.arange(1, n + 1), TYPE: np.array([lookup[e] for e in names], dtype=np.int64)}
        top = f.topology
        if top.residues:
            resid = np.zeros(n, dtype=np.int64)
            for r in top.residues:
                resid[np.asarray(r.atoms, dtype=np.int64)] = r.id if r.id is not None else 0
            props["residue"] = resid
        if f.has_velocities():
            props["velocity"] = np.array(f.velocities, dtype=float)
        if f.cell.shape == self.cf.CellShape.Infinite:
            lo, hi = (pos.min(0), pos.max(0)) if n else (np.zeros(3), np.ones(3))
            cell, origin, pbc = np.diag(np.maximum(hi - lo + 2.0, 1.0)), lo - 1.0, (False, False, False)
        else:
            cell, origin, pbc = np.array(f.cell.matrix, dtype=float).T, np.zeros(3), (True, True, True)
        fr = Frame(pos, cell, origin, pbc, props, {"Timestep": int(f.step)}, type_names={v: k for k, v in lookup.items()})
        b = np.asarray(top.bonds, dtype=np.int64).reshape(-1, 2)
        if len(b):
            fr.bonds = b
        return fr


_CHEM_READERS: dict = {}


def chem_reader(path, topology=None) -> ChemReader:
    key = (os.path.abspath(path), os.path.getmtime(path), topology)
    if key not in _CHEM_READERS:
        _CHEM_READERS.clear()
        _CHEM_READERS[key] = ChemReader(path, topology)
    return _CHEM_READERS[key]


# --------------------------------------------------------------------------
# ASE bridge
# --------------------------------------------------------------------------
def from_ase(atoms) -> Frame:
    symbols = atoms.get_chemical_symbols()
    uniq = list(dict.fromkeys(symbols))
    lookup = {e: i + 1 for i, e in enumerate(uniq)}
    props = {TYPE: np.array([lookup[s] for s in symbols], dtype=np.int64)}
    for key, arr in atoms.arrays.items():
        if key in ("positions", "numbers"):
            continue
        a = np.asarray(arr)
        if a.dtype.kind in "biuf":
            props[{"momenta": "momentum"}.get(key, key)] = a
    cell = np.array(atoms.cell)
    pbc = tuple(atoms.pbc)
    pos = atoms.positions
    origin = np.zeros(3)
    if abs(np.linalg.det(cell)) < 1e-9:
        lo, hi = pos.min(axis=0), pos.max(axis=0)
        cell = np.diag(np.maximum(hi - lo + 2.0, 1.0))
        origin = lo - 1.0
        pbc = (False, False, False)
    attrs = {k: v for k, v in atoms.info.items() if isinstance(v, (int, float, str))}
    attrs.setdefault("Timestep", 0)
    return Frame(pos, cell, origin, pbc, props, attrs, type_names={i: e for e, i in lookup.items()})


def to_ase(frame: Frame):
    from ase import Atoms
    from ase.data import chemical_symbols
    types = frame.properties.get(TYPE, np.ones(frame.n, dtype=int))
    symbols = []
    for t in types:
        name = frame.type_names.get(int(t), "")
        symbols.append(name if name in chemical_symbols else "X")
    return Atoms(symbols=symbols, positions=frame.positions - frame.origin, cell=frame.cell, pbc=frame.pbc)


# --------------------------------------------------------------------------
# Writers
# --------------------------------------------------------------------------
def _flat_columns(frame: Frame):
    names, cols = [], []
    for name, arr in frame.properties.items():
        a = np.asarray(arr)
        if a.dtype.kind not in "biuf":
            continue
        if a.ndim == 1:
            names.append(name)
            cols.append(a)
        else:
            for i in range(a.shape[1]):
                names.append(f"{name}[{i}]")
                cols.append(a[:, i])
    return names, cols


def write_lammps_dump(frame: Frame, path, append=False):
    names, cols = _flat_columns(frame)
    names = [n for n in names if n not in (ID, TYPE)]
    cols_by = dict(zip(*_flat_columns(frame)))
    ids = frame.properties.get(ID, np.arange(1, frame.n + 1))
    types = frame.properties.get(TYPE, np.ones(frame.n, dtype=int))
    a, b, c = frame.cell
    tri = abs(b[0]) + abs(c[0]) + abs(c[1]) > 1e-12 or abs(a[1]) + abs(a[2]) + abs(b[2]) > 1e-12
    lo = frame.origin
    xy, xz, yz = b[0], c[0], c[1]
    xhi, yhi, zhi = lo[0] + a[0], lo[1] + b[1], lo[2] + c[2]
    flags = " ".join("pp" if p else "ff" for p in frame.pbc)
    with open(path, "a" if append else "w") as fh:
        fh.write(f"ITEM: TIMESTEP\n{int(frame.attributes.get('Timestep', 0))}\n")
        fh.write(f"ITEM: NUMBER OF ATOMS\n{frame.n}\n")
        if tri:
            fh.write(f"ITEM: BOX BOUNDS xy xz yz {flags}\n")
            xlo_b = lo[0] + min(0, xy, xz, xy + xz)
            xhi_b = xhi + max(0, xy, xz, xy + xz)
            ylo_b = lo[1] + min(0, yz)
            yhi_b = yhi + max(0, yz)
            fh.write(f"{xlo_b} {xhi_b} {xy}\n{ylo_b} {yhi_b} {xz}\n{lo[2]} {zhi} {yz}\n")
        else:
            fh.write(f"ITEM: BOX BOUNDS {flags}\n{lo[0]} {xhi}\n{lo[1]} {yhi}\n{lo[2]} {zhi}\n")
        header = ["id", "type", "x", "y", "z"] + names
        fh.write("ITEM: ATOMS " + " ".join(h.replace(" ", "_") for h in header) + "\n")
        table = np.column_stack([ids, types, frame.positions] + [cols_by[n] for n in names])
        fmt = ["%d", "%d", "%.8g", "%.8g", "%.8g"] + ["%.8g"] * len(names)
        np.savetxt(fh, table, fmt=fmt)


def write_xyz(frame: Frame, path, append=False):
    types = frame.properties.get(TYPE, np.ones(frame.n, dtype=int))
    species = np.array([frame.type_names.get(int(t), str(int(t))) for t in types])
    names, cols = _flat_columns(frame)
    keep = [(n, c) for n, c in zip(names, cols) if n != TYPE]
    lattice = " ".join(f"{v:.8g}" for v in frame.cell.ravel())
    pbc = " ".join("T" if p else "F" for p in frame.pbc)
    props = "species:S:1:pos:R:3" + "".join(f":{n.replace(':', '_').replace('[', '_').replace(']', '')}:R:1" for n, _ in keep)
    origin = " ".join(f"{v:.8g}" for v in frame.origin)
    with open(path, "a" if append else "w") as fh:
        fh.write(f"{frame.n}\n")
        fh.write(f'Lattice="{lattice}" Origin="{origin}" Properties={props} pbc="{pbc}" '
                 f'Timestep={int(frame.attributes.get("Timestep", 0))}\n')
        for i in range(frame.n):
            extra = " ".join(f"{c[i]:.8g}" for _, c in keep)
            p = frame.positions[i]
            fh.write(f"{species[i]} {p[0]:.8g} {p[1]:.8g} {p[2]:.8g} {extra}\n")


def _savemat73(path, v):
    # MATLAB -v7.3 layout (HDF5, column-major, 512-byte MATLAB header) for data beyond the 2 GB v5 limit.
    import h5py
    with h5py.File(path, "w", userblock_size=512) as f:
        for k, a in v.items():
            a = np.asarray(a)
            if a.dtype.kind in "OU":
                continue
            d = f.create_dataset(k, data=a.T if a.ndim >= 2 else a.reshape(1, -1) if a.ndim == 1 else a)  # 1D -> MATLAB N x 1
            d.attrs["MATLAB_class"] = np.bytes_({"f": "double", "i": "int64", "u": "uint64", "b": "logical"}[a.dtype.kind])
    with open(path, "r+b") as fh:
        fh.write(b"MATLAB 7.3 MAT-file, Platform: atomview, Created by: AtomView HDF5 schema 1.00 .".ljust(116, b" ")
                 + b"\x00" * 8 + b"\x00\x02IM")


def write_mat(frames, path):
    """Frames as a .mat file in the layout read_mat understands: Position N x 3 (x F), Cell, Origin, PBC, Timestep,
    every numeric property (N x k x F), attr_<name> (one value per frame) and table_<name> ([x y], last frame)."""
    frames = [frames] if isinstance(frames, Frame) else list(frames)
    f0, F = frames[0], len(frames)
    if F > 1 and any(f.n != f0.n for f in frames): raise ValueError("Frames have different particle counts; write .xyz or one .mat per frame")
    stack = (lambda xs: np.stack(xs, axis=-1)) if F > 1 else (lambda xs: xs[0])
    cells = np.stack([f.cell for f in frames], axis=2)
    v = {"Position": stack([f.positions for f in frames]), "Cell": cells if F > 1 and np.ptp(cells, axis=2).any() else f0.cell,
         "Origin": f0.origin[None, :], "PBC": np.array(f0.pbc, dtype=float)[None, :],
         "Timestep": np.array([float(f.attributes.get("Timestep", 0) or 0) for f in frames])[None, :]}
    from .data import sanitize
    for name, a0 in f0.properties.items():
        if np.asarray(a0).dtype.kind not in "biuf" or not all(name in f.properties for f in frames):
            continue
        v[sanitize(name)] = stack([np.asarray(f.properties[name]) for f in frames])
    for k in f0.attributes:
        vals = [f.attributes.get(k) for f in frames]
        if not k.startswith("_") and all(isinstance(x, (int, float, np.integer, np.floating)) for x in vals):
            v["attr_" + sanitize(k)] = np.array(vals, dtype=float)[None, :]
    for k, t in frames[-1].tables.items():
        v[("table_" + sanitize(k))[:63]] = np.column_stack([t["x"], t["y"]])
    if f0.type_names:
        v["type_names"] = np.array([f0.type_names.get(t, "") for t in range(1, max(f0.type_names) + 1)], dtype=object)
    if sum(np.asarray(a).nbytes for a in v.values() if np.asarray(a).dtype.kind != "O") > 1.5 * 2**30:
        _savemat73(path, v)
    else:
        from scipy.io import savemat
        savemat(path, v, do_compression=True)


def write(frame: Frame, path, append=False):
    p = str(path).lower()
    if p.endswith(".mat"):
        write_mat(frame, path)
    elif p.endswith((".xyz", ".extxyz")):
        write_xyz(frame, path, append)
    elif p.endswith((".dump", ".lammpstrj")):
        write_lammps_dump(frame, path, append)
    else:
        from ase.io import write as ase_write
        ase_write(path, to_ase(frame))
