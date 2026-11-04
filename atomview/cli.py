"""Command line:

    atomview [file]                                   open the GUI
    atomview render FILE [-s session.atomview] -o fig.png [--frame 5] [--size 2400x1800] [--view iso]
    atomview animate FILE [-s session.atomview] -o movie.mp4 [--fps 15] [--every 2]
    atomview run FILE script.py [-o out.xyz]          apply a Python script modifier, write the result
    atomview info FILE                                print a summary
    atomview plugins                                  list modifiers, loaded plugins and plugin errors
    atomview msd FILE [-s session.atomview] [--every 1] [--dt 1.0] [-o msd.csv]   MSD(t) and diffusion coefficient
    atomview dynamics FILE [--k 7.0] [--every 1] [--dt 1.0] -o dyn.csv   MSD, alpha2, Fs(k,t), VACF in one CSV
    atomview convert FILE OUT.mat|.xyz|.dump [-m MOD] [--frames 0:100:5]   data or pipeline results to a file
    atomview script SCRIPT.py [args]                  run a Python script with atomview available (also in builds)
    Any command that reads a file accepts -m "Modifier(key=value)" (repeatable), e.g.
      -m "PolyhedralTemplateMatching()" -m "ColorCoding(property='c_pe', colormap='plasma')"
    atomview batch FOLDER [-o OUTDIR] [-s session.atomview] [--pattern "*.xyz,*.mat"] [--pick name|mtime]
          [--files last|all] [--frames last|all|every5|::5|-10:|0,50] [-j N]
                                                      images of the snapshots in every sub-folder
"""
from __future__ import annotations

import argparse
import os
import sys


def _pipeline(args):
    from . import modifiers  # noqa: F401
    from .io import load
    from .pipeline import Pipeline
    from .render import RenderSettings
    if getattr(args, "session", None):
        p = Pipeline.load(args.session, source=load(args.file, getattr(args, "topology", None)) if args.file else None)
        settings = RenderSettings.from_dict(p.session_extra.get("render"))
        camera = p.session_extra.get("camera")
    else:
        p = Pipeline(load(args.file, getattr(args, "topology", None)))
        settings, camera = RenderSettings(), None
    from .pipeline import parse_modifier
    for spec in getattr(args, "modifier", None) or []:
        p.add(parse_modifier(spec))
    return p, settings, camera


def _size(s):
    w, h = s.lower().split("x")
    return int(w), int(h)


def main(argv=None):
    # short error message instead of a traceback; ATOMVIEW_DEBUG=1 keeps the traceback
    if os.environ.get("ATOMVIEW_DEBUG"): return _main(argv)
    try: return _main(argv)
    except Exception as e:
        print(f"error: {type(e).__name__}: {e}", file=sys.stderr); return 1


def _main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    from .plugins import ERRORS, LOADED, load_plugins
    load_plugins()
    for name, err in ERRORS:
        print(f"warning: plugin {name} failed to load:\n{err}", file=sys.stderr)
    if argv[:1] == ["plugins"]:
        from .pipeline import REGISTRY
        print("loaded plugins:", *(LOADED or ["(none)"]), sep="\n  ")
        print("modifiers:", *sorted(f"{c.category:12s} {c.title}  [{n}]" for n, c in REGISTRY.items()), sep="\n  ")
        return 1 if ERRORS else 0
    commands = {"render", "animate", "run", "info", "batch", "msd", "dynamics", "convert", "script"}
    if not argv or argv[0] not in commands:
        if argv and argv[0] in ("-h", "--help"):
            print(__doc__)
            return 0
        from .gui.main_window import run
        return run(argv[0] if argv else None)

    ap = argparse.ArgumentParser(prog="atomview")
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("render", "animate"):
        s = sub.add_parser(name)
        s.add_argument("file", nargs="?")
        s.add_argument("-s", "--session")
        s.add_argument("-m", "--modifier", action="append", help="e.g. \"ColorCoding(property='c_pe')\" (repeatable)")
        s.add_argument("-t", "--topology", help="topology for XTC/TRR/DCD/NetCDF (.gro, .pdb, .psf, .tpr)")
        s.add_argument("-o", "--output", required=True)
        s.add_argument("--size", type=_size, default=(2400, 1800) if name == "render" else (1280, 960))
        s.add_argument("--view", default="iso", choices=["iso", "+x", "-x", "+y", "-y", "+z", "-z"])
        s.add_argument("--roll", type=float, default=None, help="rotate the view in the image plane (degrees)")
        s.add_argument("--ao", action="store_true", help="ambient occlusion")
        s.add_argument("--background", default=None, help="white, black or r,g,b")
        s.add_argument("--overlay", default=None, help='per-frame text, e.g. "t = {Timestep*0.002:.1f} ps"')
        if name == "render":
            s.add_argument("--frame", type=int, default=0)
            s.add_argument("--transparent", action="store_true")
        else:
            s.add_argument("--fps", type=int, default=10)
            s.add_argument("--every", type=int, default=1)
    r = sub.add_parser("run")
    r.add_argument("file")
    r.add_argument("script")
    r.add_argument("-o", "--output")
    r.add_argument("--frame", type=int, default=None, help="only this frame (default: all)")
    i = sub.add_parser("info")
    i.add_argument("file")
    i.add_argument("-t", "--topology")
    b = sub.add_parser("batch", help="render the last file of every folder under FOLDER")
    b.add_argument("folder")
    b.add_argument("-o", "--output", help="put all images here (default: next to each file)")
    b.add_argument("-s", "--session", help="session JSON whose modifiers and render settings are applied")
    b.add_argument("-m", "--modifier", action="append", help="e.g. \"ColorCoding(property='c_pe')\" (repeatable)")
    b.add_argument("--pattern", default="*.xyz,*.mat", help='comma-separated, e.g. "*.xyz,*.mat,*.dump"')
    b.add_argument("--pick", default="name", choices=["name", "mtime"],
                   help="last by natural file name (default) or by modification time")
    b.add_argument("--frame", type=int, default=-1, help="frame within the file; -1 = last (default)")
    b.add_argument("--frames", default=None,
                   help='which snapshots per folder: last (default), first, all, every5, "::5", "--frames=-10:", "0,50,100"')
    b.add_argument("--files", default="last", choices=["last", "all"],
                   help="snapshots from the last matching file (default) or from all matching files in order")
    b.add_argument("-j", "--jobs", type=int, default=1, help="parallel processes; 0 = one per core")
    b.add_argument("--size", type=_size, default=(2000, 1500))
    b.add_argument("--view", default="iso", choices=["iso", "+x", "-x", "+y", "-y", "+z", "-z"])
    b.add_argument("--roll", type=float, default=None, help="rotate the view in the image plane (degrees)")
    b.add_argument("--format", default="png", choices=["png", "jpg", "tif", "svg", "pdf"])
    b.add_argument("--no-recursive", action="store_true", help="only FOLDER itself")
    b.add_argument("--skip-existing", action="store_true", help="skip folders whose image is newer than the file")
    b.add_argument("--transparent", action="store_true")
    b.add_argument("--ao", action="store_true", help="ambient occlusion")
    b.add_argument("--overlay", default=None, help='text on each image, e.g. "step {Timestep}  N = {N}"')
    b.add_argument("--session-camera", action="store_true", help="use the camera saved in the session")
    q = sub.add_parser("msd", help="mean squared displacement over the trajectory")
    q.add_argument("file")
    q.add_argument("-s", "--session")
    q.add_argument("-m", "--modifier", action="append", help="e.g. \"ColorCoding(property='c_pe')\" (repeatable)")
    q.add_argument("-t", "--topology", help="topology for XTC/TRR/DCD/NetCDF")
    q.add_argument("-o", "--output", help="CSV file (time, MSD)")
    q.add_argument("--every", type=int, default=1)
    q.add_argument("--dt", type=float, default=None, help="time per frame (default: from Timestep)")
    q.add_argument("--selected", action="store_true", help="only particles selected in the first frame")
    y = sub.add_parser("dynamics", help="MSD, non-Gaussian parameter, Fs(k,t) and VACF over the trajectory")
    y.add_argument("file")
    y.add_argument("-s", "--session")
    y.add_argument("-m", "--modifier", action="append", help="e.g. \"ColorCoding(property='c_pe')\" (repeatable)")
    y.add_argument("-t", "--topology")
    y.add_argument("-o", "--output", required=True, help="CSV file")
    y.add_argument("--k", type=float, default=None, help="wave number for Fs(k,t)")
    y.add_argument("--every", type=int, default=1)
    y.add_argument("--dt", type=float, default=None)
    y.add_argument("--selected", action="store_true")
    c = sub.add_parser("convert", help="write a file (after any modifiers) as .mat, .xyz or .dump")
    c.add_argument("file")
    c.add_argument("output")
    c.add_argument("-s", "--session")
    c.add_argument("-t", "--topology")
    c.add_argument("-m", "--modifier", action="append", help="e.g. \"ColorCoding(property='c_pe')\" (repeatable)")
    c.add_argument("--frames", default=None, help="start:stop[:step] (default: all)")
    k = sub.add_parser("script", help="run a Python script with atomview importable")
    k.add_argument("script")
    k.add_argument("script_args", nargs=argparse.REMAINDER)
    args = ap.parse_args(argv)

    if args.cmd == "script":
        import runpy
        sys.argv = [args.script, *args.script_args]
        sys.path.insert(0, os.path.dirname(os.path.abspath(args.script)))
        runpy.run_path(args.script, run_name="__main__")
        return 0

    if args.cmd == "convert":
        from .io import write, write_mat
        p, _, _ = _pipeline(args)
        sl = slice(*[int(x) if x else None for x in args.frames.split(":")]) if args.frames else slice(None)
        frames = list(range(p.n_frames))[sl]
        for name, err in p.errors():
            print(f"warning: {name}: {err}", file=sys.stderr)
        if args.output.lower().endswith(".mat"):
            write_mat([p.compute(f) for f in frames], args.output)
        else:
            for n, f in enumerate(frames):
                write(p.compute(f), args.output, append=n > 0)
        for name, err in p.errors():
            print(f"warning: {name}: {err}", file=sys.stderr)
        print(args.output)
        return 0

    if args.cmd == "dynamics":
        import numpy as np
        from .trajectory import dynamics
        p, _, _ = _pipeline(args)
        res = dynamics(p, range(0, p.n_frames, args.every), args.dt, args.k, args.selected)
        names = [n for n in res if n != "VACF(t)"]
        F = len(res["MSD(t)"]["x"])
        vacf = np.full(F, np.nan)
        vacf[:len(res["VACF(t)"]["y"])] = res["VACF(t)"]["y"]  # finite differences give one point fewer
        cols = [res["MSD(t)"]["x"]] + [res[n]["y"] for n in names] + [vacf]
        header = ",".join(["time"] + [n.replace(",", ";") for n in names] + ["VACF(t)"])
        np.savetxt(args.output, np.column_stack(cols), delimiter=",", header=header, comments="")
        r = res["MSD(t)"]
        print(f"{r['frames']} frames, {r['particles']} particles: D = {r['D']:.6g}; VACF from {res['VACF(t)']['source']}")
        print(args.output)
        return 0

    if args.cmd == "msd":
        import numpy as np
        from .trajectory import msd
        p, _, _ = _pipeline(args)
        r = msd(p, range(0, p.n_frames, args.every), args.dt, args.selected)
        print(f"{r['frames']} frames, {r['particles']} particles, {r['dim']}D: D = {r['D']:.6g} (linear fit over the second half)"
              + (f"; {r['dropped']} particles missing in some frame were left out" if r['dropped'] else ""))
        if args.output:
            np.savetxt(args.output, np.column_stack([r["x"], r["y"]]), delimiter=",", header="time,msd", comments="")
            print(args.output)
        return 0

    if args.cmd == "batch":
        from . import batch
        from .render import RenderSettings
        settings = None
        if args.ao or args.overlay is not None or args.roll is not None:
            import json
            saved = {}
            if args.session:
                with open(args.session) as fh:
                    saved = json.load(fh).get("render", {})
            settings = RenderSettings.from_dict(saved)
            settings.ambient_occlusion = settings.ambient_occlusion or args.ao
            if args.overlay is not None:
                settings.overlay = args.overlay
            if args.roll is not None:
                settings.view_roll = args.roll
        from .pipeline import parse_modifier
        extra = [parse_modifier(m).to_dict() for m in args.modifier or []]
        results = batch.run(args.folder, patterns=[p.strip() for p in args.pattern.split(",") if p.strip()],
                            out_dir=args.output, session=args.session, modifiers=extra, settings=settings, pick=args.pick,
                            frame=args.frame, size=args.size, view=args.view, fmt=args.format,
                            recursive=not args.no_recursive, skip_existing=args.skip_existing,
                            transparent=args.transparent, use_session_camera=args.session_camera,
                            frames=args.frames, files=args.files, jobs=args.jobs)
        failed = [r for r in results if r["status"].startswith("error")]
        skipped = sum(r["status"].startswith("skipped") for r in results)
        print(f"{len(results) - len(failed) - skipped} image(s) written, {skipped} up to date, {len(failed)} failed")
        return 1 if failed else 0

    if args.cmd == "info":
        from .io import load
        t = load(args.file, args.topology)
        f = t.get(0)
        print(f"{args.file}: format {t.format}, {t.n_frames} frame(s)")
        print(f.summary())
        print("attributes:", {k: v for k, v in f.attributes.items() if not k.startswith("_")})
        return 0

    if args.cmd == "run":
        from .io import load, write
        from .modifiers import MatlabScript, PythonScript
        from .pipeline import Pipeline
        p = Pipeline(load(args.file))
        p.strict = True
        p.add(MatlabScript(file=args.script) if args.script.endswith(".m") else PythonScript(file=args.script))
        frames = [args.frame] if args.frame is not None else range(p.n_frames)
        for k, fr in enumerate(frames):
            out = p.compute(fr)
            if p.modifiers[0].log:
                print(p.modifiers[0].log)
            if args.output:
                write(out, args.output, append=k > 0)
        return 0

    p, settings, camera = _pipeline(args)
    if args.ao:
        settings.ambient_occlusion = True
    if args.overlay is not None:
        settings.overlay = args.overlay
    if args.roll is not None:
        settings.view_roll = args.roll
    if args.background:
        named = {"white": (1, 1, 1), "black": (0, 0, 0)}
        settings.background = named.get(args.background) or tuple(float(v) for v in args.background.split(","))
    if p.errors():
        for name, err in p.errors():
            print(f"warning: {name}: {err}", file=sys.stderr)
    if args.cmd == "render":
        from .render import render_image
        render_image(p, args.output, frame=args.frame, size=args.size, settings=settings, camera=camera,
                     transparent=args.transparent, view=args.view)
    else:
        from .render import render_animation

        def prog(i, n):
            print(f"\rframe {i}/{n}", end="", file=sys.stderr)

        render_animation(p, args.output, frames=range(0, p.n_frames, args.every), size=args.size, fps=args.fps,
                         settings=settings, camera=camera, view=args.view, progress=prog)
        print(file=sys.stderr)
    for name, err in p.errors():
        print(f"warning: {name}: {err}", file=sys.stderr)
    print(args.output)
    return 0


if __name__ == "__main__":
    sys.exit(main())
