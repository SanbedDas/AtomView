"""Pipelines: a data source followed by an ordered list of modifiers."""
from __future__ import annotations

import glob
import json
import os
import time
from dataclasses import dataclass, field
from typing import Any

from .data import Frame

REGISTRY: dict[str, type["Modifier"]] = {}


class Cancelled(Exception):
    pass


@dataclass
class Param:
    name: str
    kind: str  # float, int, bool, choice, str, expr, code, property, color, file, vector
    default: Any
    label: str = ""
    choices: tuple = ()
    min: float | None = None
    max: float | None = None
    step: float | None = None
    help: str = ""

    def __post_init__(self):
        if not self.label:
            self.label = self.name.replace("_", " ").capitalize()


class Modifier:
    """Base class. Subclasses declare `title`, `category`, `params` and implement apply()."""

    title = "Modifier"
    category = "Other"
    description = ""
    params: list[Param] = []

    def __init_subclass__(cls, **kw):
        super().__init_subclass__(**kw)
        if cls.__dict__.get("title") and not cls.__name__.startswith("_"):
            REGISTRY[cls.__name__] = cls

    def __init__(self, **values):
        self.enabled = True
        self.values = {p.name: p.default for p in self.params}
        unknown = set(values) - set(self.values)
        if unknown:
            raise TypeError(f"{type(self).__name__} has no parameter(s) {', '.join(sorted(unknown))}")
        self.values.update(values)
        self.status = ""  # message shown under the modifier after evaluation
        self.error = None
        self.log = ""

    def __getattr__(self, name):
        values = self.__dict__.get("values")
        if values is not None and name in values:
            return values[name]
        raise AttributeError(name)

    def set(self, **values):
        self.values.update(values)
        return self

    def label(self) -> str:
        return self.title

    def apply(self, frame: Frame, ctx: "Context") -> Frame:  # pragma: no cover - abstract
        raise NotImplementedError

    def to_dict(self):
        return {"type": type(self).__name__, "enabled": self.enabled, "values": dict(self.values)}

    @staticmethod
    def from_dict(d):
        cls = REGISTRY[d["type"]]
        m = cls()
        m.values.update({k: v for k, v in d.get("values", {}).items() if k in m.values})
        m.enabled = d.get("enabled", True)
        return m


@dataclass
class Context:
    frame_index: int
    pipeline: "Pipeline"
    upto: int = 0  # index of the modifier being applied

    def cancelled(self) -> bool:
        return self.pipeline.cancel.is_set()

    def check(self):
        if self.pipeline.cancel.is_set(): raise Cancelled()

    def upstream(self, frame_index: int) -> Frame:
        """The pipeline output just before the current modifier, at another frame."""
        return self.pipeline.compute(frame_index, stop=self.upto)


@dataclass
class Pipeline:
    source: Any  # anything with n_frames and get(i) -> Frame (e.g. io.Trajectory)
    modifiers: list = field(default_factory=list)

    def __post_init__(self):
        import threading
        self._cache = {}
        self.cancel = threading.Event()  # set from another thread to stop a running compute()
        self.progress = None  # optional callback(message) called before each modifier
        self.revision = 0
        self.last_compute_time = 0.0

    @property
    def n_frames(self):
        return self.source.n_frames

    def add(self, modifier) -> "Modifier":
        if callable(modifier) and not isinstance(modifier, Modifier):
            from .modifiers import Function
            modifier = Function(modifier)
        self.modifiers.append(modifier)
        self.invalidate()
        return modifier

    def invalidate(self):
        self._cache.clear()
        self.revision += 1

    def compute(self, frame_index: int = 0, stop: int | None = None) -> Frame:
        stop = len(self.modifiers) if stop is None else stop
        key = (frame_index, stop)
        if key in self._cache:
            return self._cache[key]
        t0 = time.perf_counter()
        frame = self.source.get(frame_index).copy()
        for k, m in enumerate(self.modifiers[:stop]):
            if not m.enabled:
                m.status = "disabled"
                continue
            if self.cancel.is_set(): raise Cancelled()
            if self.progress: self.progress(f"{m.label()} ({k + 1}/{stop})")
            try:
                m.error = None
                frame = m.apply(frame, Context(frame_index, self, k)) or frame
            except Cancelled:
                raise
            except Exception as e:  # keep going so the viewer still shows something
                m.error = f"{type(e).__name__}: {e}"
                m.status = m.error
                if getattr(self, "strict", False):
                    raise
        if self.cancel.is_set(): raise Cancelled()
        if stop == len(self.modifiers):
            self.last_compute_time = time.perf_counter() - t0
        if len(self._cache) > 16:
            self._cache.pop(next(iter(self._cache)))
        self._cache[key] = frame
        return frame

    def errors(self):
        return [(m.label(), m.error) for m in self.modifiers if m.error]

    # ---- sessions ------------------------------------------------------------
    def to_dict(self):
        return {"source": getattr(self.source, "path", None), "topology": getattr(self.source, "topology", None),
                "modifiers": [m.to_dict() for m in self.modifiers]}

    def save(self, path, extra=None):
        d = {"atomview_session": SESSION_VERSION, **self.to_dict()}
        for k in ("source", "topology"):  # also relative to the session, so a moved project folder still opens
            d[k + "_relative"] = _relpath(d[k], os.path.dirname(os.path.abspath(path))) if d[k] else None
        if extra:
            d.update(extra)
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(d, fh, indent=2, default=_json_default)

    @staticmethod
    def load(path, source=None):
        from . import modifiers  # noqa: F401  (registers modifier classes)
        from .plugins import load_plugins
        load_plugins()
        from .io import load as load_traj
        with open(path, encoding="utf-8") as fh:
            d = json.load(fh)
        if source is None:
            source = load_traj(session_file(d, path, "source"), session_file(d, path, "topology"))
        p = Pipeline(source, [Modifier.from_dict(m) for m in d.get("modifiers", [])])
        p.session_extra = {k: v for k, v in d.items() if k not in ("atomview_session", "source", "topology", "source_relative", "topology_relative", "modifiers")}
        return p


SESSION_VERSION = 1


def _relpath(f, base):
    try:
        return os.path.relpath(f, base)
    except ValueError:  # another drive on Windows
        return None


def _found(f):
    return os.path.exists(f) or (any(c in f for c in "*?[") and bool(glob.glob(f)))


def session_file(d, path, key):
    """A file named in a session: the saved path, else the same file relative to the session file."""
    f, rel = d.get(key), d.get(key + "_relative")
    if not f:
        return None
    if _found(f):
        return f
    if rel:
        g = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(path)), rel))
        if _found(g):
            return g
    raise FileNotFoundError(f"The {'data' if key == 'source' else key} file of this session was not found: {f}")


def is_session(path):
    if path.lower().endswith(".atomview"):
        return True
    if not path.lower().endswith(".json"):
        return False
    try:
        with open(path, encoding="utf-8") as fh:
            d = json.load(fh)
        return isinstance(d, dict) and "modifiers" in d
    except (OSError, ValueError):
        return False


def _json_default(o):
    import numpy as np
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    return str(o)


def parse_modifier(spec: str) -> "Modifier":
    """A modifier from text such as "ColorCoding(property='c_pe', colormap='plasma')": a registered class name
    with keyword arguments that are Python literals. Nothing is executed."""
    import ast
    from . import modifiers  # noqa: F401
    if "(" not in spec:  # bare class name or title, e.g. "Polyhedral template matching"
        cls = REGISTRY.get(spec.strip()) or {c.title.lower(): c for c in REGISTRY.values()}.get(spec.strip().lower())
        if cls is None: raise ValueError(f"Unknown modifier {spec!r}. Known: {', '.join(sorted(REGISTRY))}")
        return cls()
    tree = ast.parse(spec.strip(), mode="eval").body
    if isinstance(tree, ast.Name):
        tree = ast.Call(func=tree, args=[], keywords=[])
    if not isinstance(tree, ast.Call) or not isinstance(tree.func, ast.Name) or tree.args:
        raise ValueError(f"Expected Name(key=value, ...), got {spec!r}")
    name = tree.func.id
    cls = REGISTRY.get(name) or {c.title.lower(): c for c in REGISTRY.values()}.get(name.lower())
    if cls is None:
        raise ValueError(f"Unknown modifier {name!r}. Known: {', '.join(sorted(REGISTRY))}")
    return cls(**{k.arg: ast.literal_eval(k.value) for k in tree.keywords})
