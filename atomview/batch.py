"""Batch rendering: walk a folder tree and render images of the snapshots in each folder.

Per folder the snapshots are the frames of the last matching file (files="last") or of all matching files in
order (files="all"). "Last" means the highest file name in natural order (run_9.xyz < run_10.xyz) or, with
pick="mtime", the most recently modified file. `frames` picks among the snapshots: "last" (default), "first",
"all", a slice such as "::5" or "-10:", "every5", an index, or a comma list "0,10,20".
"""
from __future__ import annotations

import csv
import fnmatch
import os
import time
import traceback

from .io import _natural_key, load
from .pipeline import Modifier, Pipeline

DEFAULT_PATTERNS = ("*.xyz", "*.mat")


def matching_files(folder: str, patterns=DEFAULT_PATTERNS, pick: str = "name") -> list:
    try:
        names = [n for n in os.listdir(folder) if os.path.isfile(os.path.join(folder, n))]
    except OSError:
        return []
    pats = [p.lower() for p in patterns]
    matches = [n for n in names if not n.startswith(".") and any(fnmatch.fnmatch(n.lower(), p) for p in pats)]
    key = (lambda n: (os.path.getmtime(os.path.join(folder, n)), _natural_key(n))) if pick == "mtime" else _natural_key
    return [os.path.join(folder, n) for n in sorted(matches, key=key)]


def latest_file(folder: str, patterns=DEFAULT_PATTERNS, pick: str = "name") -> str | None:
    fs = matching_files(folder, patterns, pick)
    return fs[-1] if fs else None


def single(spec) -> bool:
    s = str(spec).strip().lower()
    return s in ("", "last", "first") or s.lstrip("-").isdigit()


def frame_spec(spec, n: int) -> list:
    """Indices into a sequence of n snapshots."""
    s = str(spec).strip().lower().replace(" ", "")
    if s in ("", "last"): return [n - 1]
    if s == "first": return [0]
    if s == "all": return list(range(n))
    if s.startswith("every"): return list(range(0, n, int(s[5:].lstrip(":=") or 1)))
    if ":" in s: return list(range(n))[slice(*[int(x) if x else None for x in s.split(":")])]
    if s.lstrip("-").isdigit():
        k = int(s); return [max(0, min(k if k >= 0 else n + k, n - 1))]
    return sorted({int(x) % n for x in s.split(",") if x and -n <= int(x) < n})


def folders(root: str, recursive: bool = True, skip=()):
    """The root and (optionally) every sub-folder, sorted, skipping hidden ones and `skip`."""
    root = os.path.abspath(root)
    skip = {os.path.abspath(s) for s in skip if s}
    if not recursive:
        yield root
        return
    for dirpath, dirnames, _ in os.walk(root):
        dirnames[:] = sorted((d for d in dirnames if not d.startswith(".")
                              and os.path.abspath(os.path.join(dirpath, d)) not in skip), key=_natural_key)
        yield dirpath


def image_path(file: str, root: str, out_dir: str | None, fmt: str = "png") -> str:
    stem = os.path.splitext(os.path.basename(file))[0]
    if not out_dir:
        return os.path.join(os.path.dirname(file), f"{stem}.{fmt}")
    rel = os.path.relpath(os.path.dirname(file), os.path.abspath(root))
    prefix = "" if rel == "." else rel.replace(os.sep, "__").replace("/", "__") + "__"
    return os.path.join(out_dir, f"{prefix}{stem}.{fmt}")


def _folder(folder, root, c, report=None):
    """Render the selected snapshots of one folder; returns one record per image."""
    from .render import render_image
    fs = matching_files(folder, c["patterns"], c["pick"])
    fs = fs[-1:] if c["files"] == "last" else fs
    many = not single(c["frames"])
    recs, seq, trajs = [], [], {}
    spec = str(c["frames"]).strip().lower()
    if not many and c["skip_existing"] and (c["files"] == "last" or spec in ("", "last", "first")):
        f = fs[0] if spec == "first" else fs[-1]  # one image named after its file: check before reading anything
        out = image_path(f, root, c["out_dir"], c["fmt"])
        if os.path.exists(out) and os.path.getmtime(out) >= os.path.getmtime(f):
            rec = dict(folder=folder, file=f, image=out, frame="", particles="", timestep="",
                       status="skipped (up to date)", seconds="0.00")
            report and report(rec)
            return [rec]
    for f in fs:
        try:
            trajs[f] = load(f)
            seq += [(f, k) for k in range(trajs[f].n_frames)]
        except Exception as e:
            rec = dict(folder=folder, file=f, image="", frame="", particles="", timestep="",
                       status=f"error: {type(e).__name__}: {e}", seconds="")
            recs.append(rec); report and report(rec)
    pipes = {}
    for f, k in [seq[i] for i in frame_spec(c["frames"], len(seq))] if seq else []:
        out = image_path(f, root, c["out_dir"], c["fmt"])
        if many and trajs[f].n_frames > 1:
            out = f"{os.path.splitext(out)[0]}_{k:05d}.{c['fmt']}"
        rec = dict(folder=folder, file=f, image=out, frame=k, particles="", timestep="", status="", seconds="")
        t0 = time.perf_counter()
        try:
            if c["skip_existing"] and os.path.exists(out) and os.path.getmtime(out) >= os.path.getmtime(f):
                rec["status"] = "skipped (up to date)"
            else:
                if f not in pipes:
                    pipes.clear()
                    pipes[f] = Pipeline(trajs[f], [Modifier.from_dict(m) for m in c["modifiers"]])
                pipe = pipes[f]
                data = pipe.compute(k)
                render_image(data, out, size=c["size"], settings=c["settings"], camera=c["camera"],
                             transparent=c["transparent"], view=c["view"])
                errs = pipe.errors()
                rec.update(particles=data.n, timestep=data.attributes.get("Timestep", ""),
                           status="ok" if not errs else "ok with warnings: " + "; ".join(f"{a}: {b}" for a, b in errs))
        except Exception as e:
            rec["status"] = f"error: {type(e).__name__}: {e}"
            if os.environ.get("ATOMVIEW_DEBUG"):
                traceback.print_exc()
        rec["seconds"] = f"{time.perf_counter() - t0:.2f}"
        recs.append(rec); report and report(rec)
    return recs


def _worker_init():
    from . import modifiers  # noqa: F401
    from .plugins import load_plugins
    load_plugins()
    try:
        import freud
        freud.parallel.set_num_threads(1)  # the processes already use every core
    except Exception:
        pass


def run(root, patterns=DEFAULT_PATTERNS, out_dir=None, session=None, modifiers=None, settings=None,
        pick="name", frame=-1, size=(2000, 1500), view="iso", fmt="png", recursive=True, skip_existing=False,
        transparent=False, use_session_camera=False, progress=print, cancel=None, frames=None, files="last", jobs=1):
    """Render images for every folder. Returns a list of result dicts (also written to batch_summary.csv).

    modifiers: list of Modifier.to_dict() entries (or a session file via `session`) applied to every file.
    frames/files: which snapshots (see the module doc); `frame` is the older single-index form.
    jobs: number of processes (0 = one per core). progress(message) is called per image; cancel() True stops.
    """
    from .render import RenderSettings

    camera = None
    if session:
        import json
        with open(session, encoding="utf-8") as fh:
            sd = json.load(fh)
        modifiers = sd.get("modifiers", []) + list(modifiers or [])
        if settings is None:
            settings = RenderSettings.from_dict(sd.get("render"))
        if use_session_camera:
            camera = sd.get("camera")
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    c = dict(patterns=list(patterns), pick=pick, files=files, frames=frames if frames is not None else
             ("last" if frame == -1 else str(frame)), out_dir=out_dir, modifiers=list(modifiers or []),
             settings=settings or RenderSettings(), camera=camera, size=tuple(size), view=view, fmt=fmt,
             skip_existing=skip_existing, transparent=transparent)
    todo = [d for d in folders(root, recursive, skip=[out_dir]) if matching_files(d, patterns)]
    if not todo:
        progress(f"No files matching {', '.join(patterns)} under {root}")
    done, out = [0], [None] * len(todo)

    def report(rec):
        done[0] += 1
        progress(f"[{done[0]}] {os.path.relpath(rec['file'], root)} -> {rec['image']}  ({rec['status']})")

    jobs = min(jobs or os.cpu_count() or 1, 61 if os.name == "nt" else 1024)  # Windows allows at most 61 workers
    if jobs == 1 or len(todo) < 2:
        for i, d in enumerate(todo):
            if cancel and cancel():
                break
            out[i] = _folder(d, root, c, report)
    else:
        import multiprocessing as mp
        from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
        with ProcessPoolExecutor(min(jobs, len(todo)), mp.get_context("spawn"), _worker_init) as ex:
            futs = {ex.submit(_folder, d, root, c): i for i, d in enumerate(todo)}
            pending = set(futs)
            while pending:
                fin, pending = wait(pending, timeout=0.5, return_when=FIRST_COMPLETED)
                for fu in fin:
                    try:
                        out[futs[fu]] = fu.result()
                    except Exception as e:  # a worker that died (e.g. out of memory) loses only its folder
                        d = todo[futs[fu]]
                        out[futs[fu]] = [dict(folder=d, file=d, image="", frame="", particles="", timestep="",
                                              status=f"error: {type(e).__name__}: {e}", seconds="")]
                    for rec in out[futs[fu]]:
                        report(rec)
                if cancel and cancel():
                    ex.shutdown(cancel_futures=True)
                    break
    results = [r for rs in out if rs for r in rs]
    summary = os.path.join(out_dir or os.path.abspath(root), "batch_summary.csv")
    if results:
        with open(summary, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(results[0]))
            w.writeheader()
            w.writerows(results)
    return results
