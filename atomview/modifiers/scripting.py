"""User scripts in the pipeline: Python functions and MATLAB / GNU Octave .m files."""
from __future__ import annotations

import contextlib
import sys
import threading
import glob
import io
import os
import shutil
import subprocess
import tempfile
import time
import traceback

import numpy as np

from ..data import SELECTION, Frame, sanitize
from ..pipeline import Cancelled, Modifier, Param

PYTHON_TEMPLATE = '''\
# Runs for every frame. `data` is the particle Frame; change it in place.
#   data.positions (N x 3)   data['type']   data['c_pe']   data.n
#   data.cell (rows = a, b, c)   data.origin   data.attributes   data.tables
# Return nothing, or return a new Frame (e.g. data.subset(mask)).
# numpy is available as np; import anything else you have installed.

def modify(frame, data):
    z = data.positions[:, 2]
    data['height'] = z - z.min()                  # new property
    data['selection'] = z > np.percentile(z, 90)  # select the top 10 %
    data.attributes['MeanHeight'] = float(data['height'].mean())
'''

MATLAB_TEMPLATE = '''\
% Runs once per frame in MATLAB or GNU Octave. Variables you get:
%   Position (N x 3), type, id and every other property by name (N x 1 or N x k)
%   Cell (3 x 3, rows = cell vectors), Origin (1 x 3), Timestep, Frame, N
% Any variable you leave with N rows becomes a particle property.
% Scalars become global attributes. Set DeleteMask (N x 1 logical) to remove particles.

c = mean(Position, 1);
r = sqrt(sum((Position - c).^2, 2));
distance_from_center = r;
selection = r < 0.5 * min(sqrt(sum(Cell.^2, 2)));
MeanRadius = mean(r);
'''


class _ThreadStdout:
    # sys.stdout replacement that sends writes from threads with a registered buffer there, and everything else
    # to the real stdout; redirect_stdout would capture every thread's output while a script runs.
    def __init__(self, real): self.real, self.bufs = real, {}

    def write(self, s): return self.bufs.get(threading.get_ident(), self.real).write(s)

    def flush(self): self.bufs.get(threading.get_ident(), self.real).flush()

    def __getattr__(self, name): return getattr(self.real, name)


_lock = threading.Lock()


@contextlib.contextmanager
def capture_stdout(buf):
    with _lock:
        if not isinstance(sys.stdout, _ThreadStdout):
            sys.stdout = _ThreadStdout(sys.stdout)
        sys.stdout.bufs[threading.get_ident()] = buf
    try:
        yield buf
    finally:
        sys.stdout.bufs.pop(threading.get_ident(), None)


class PythonScript(Modifier):
    title = "Python script"
    category = "Scripting"
    description = "Your own Python code: def modify(frame, data). Full numpy/scipy/freud/ASE access."
    params = [Param("code", "code", PYTHON_TEMPLATE),
              Param("file", "file", "", help="Optional .py file; overrides the code above")]

    def __init__(self, **kw):
        super().__init__(**kw)
        self._compiled_for = None
        self._func = None

    def label(self):
        return f"Python script: {os.path.basename(self.file)}" if self.file else "Python script"

    def _source(self):
        if self.file:
            with open(os.path.expanduser(self.file)) as fh:
                return fh.read()
        return self.code

    def _load(self):
        src = self._source()
        if src != self._compiled_for:
            ns = {"np": np, "__name__": "__atomview_script__", "Frame": Frame}
            exec(compile(src, self.file or "<script>", "exec"), ns)
            func = ns.get("modify")
            if not callable(func):
                raise ValueError("The script must define a function modify(frame, data)")
            self._func, self._compiled_for = func, src
        return self._func

    def apply(self, frame, ctx):
        func = self._load()
        out = io.StringIO()
        try:
            with capture_stdout(out):
                result = func(ctx.frame_index, frame)
        except Exception:
            self.log = out.getvalue() + traceback.format_exc(limit=-3)
            raise
        self.log = out.getvalue()
        self.status = "ok" + (f" — {len(self.log.splitlines())} lines of output" if self.log else "")
        return result if isinstance(result, Frame) else frame


# --------------------------------------------------------------------------
# MATLAB / Octave
# --------------------------------------------------------------------------
_MATLAB_GLOBS = (
    "/Applications/MATLAB_R*.app/bin/matlab",                 # macOS
    "/usr/local/MATLAB/R*/bin/matlab", "/opt/MATLAB/R*/bin/matlab", os.path.expanduser("~/MATLAB/R*/bin/matlab"),
    r"C:\Program Files\MATLAB\R*\bin\matlab.exe",              # Windows
)
_OCTAVE_GLOBS = (
    "/opt/homebrew/bin/octave-cli", "/usr/local/bin/octave-cli", "/usr/bin/octave-cli",
    "/Applications/Octave*.app/Contents/Resources/usr/bin/octave-cli",
    r"C:\Program Files\GNU Octave\Octave-*\mingw64\bin\octave-cli.exe",
    os.path.expanduser(r"~\AppData\Local\Programs\GNU Octave\Octave-*\mingw64\bin\octave-cli.exe"),
)


def _newest(patterns):
    hits = [p for pat in patterns for p in glob.glob(pat)]
    return sorted(hits)[-1] if hits else None


def find_engines():
    """Available engines as {'matlab': path, 'octave': path} on macOS, Linux or Windows.

    Set ATOMVIEW_MATLAB / ATOMVIEW_OCTAVE to point at a specific executable.
    """
    found = {}
    m = os.environ.get("ATOMVIEW_MATLAB") or shutil.which("matlab") or _newest(_MATLAB_GLOBS)
    if m:
        found["matlab"] = m
    o = (os.environ.get("ATOMVIEW_OCTAVE") or shutil.which("octave-cli") or shutil.which("octave")
         or _newest(_OCTAVE_GLOBS))
    if o:
        found["octave"] = o
    return found


def _mq(path: str) -> str:
    """A path as a MATLAB single-quoted string literal (forward slashes work on every OS)."""
    return "'" + str(path).replace("\\", "/").replace("'", "''") + "'"


def _mat_name(name: str) -> str:
    s = sanitize(name)
    if not s[0].isalpha():
        s = "p_" + s
    return s[:63]


_SERVER_M = r"""function av_server(d)
fid = fopen(fullfile(d, 'ready'), 'w'); fclose(fid);
while true
  j = fullfile(d, 'job.txt');
  if exist(j, 'file') == 2
    s = strtrim(fileread(j)); delete(j);
    if strcmp(s, 'quit'), break; end
    L = strsplit(s, char(10));
    err = '';
    try
      av_job(L{1}, L{2});
    catch e
      err = [e.identifier ': ' e.message];
    end
    fid = fopen(fullfile(L{1}, 'done.tmp'), 'w'); fprintf(fid, '%s', err); fclose(fid);
    movefile(fullfile(L{1}, 'done.tmp'), fullfile(L{1}, 'done.txt'));
  end
  pause(0.02);
end
end
"""

_JOB_M = r"""function av_job(av_dir, av_script)
cd(av_dir);
diary(fullfile(av_dir, 'log.txt'));
load('in.mat');
[~, av_name] = fileparts(av_script);
clear(av_name);
rehash;
run(av_script);
diary off;
cd(av_dir);
clear ans av_dir av_script av_name;
SAVE
end
"""


class _Session:
    # One long-running MATLAB/Octave process executing jobs from a folder: in.mat + script -> out.mat.
    def __init__(self, kind, exe):
        self.kind, self.exe, self.n = kind, exe, 0
        self.dir = tempfile.mkdtemp(prefix="atomview_session_")
        with open(os.path.join(self.dir, "av_server.m"), "w") as fh: fh.write(_SERVER_M)
        save = "save('-v7', 'out.mat');" if kind == "octave" else "save('out.mat', '-v7');"
        with open(os.path.join(self.dir, "av_job.m"), "w") as fh: fh.write(_JOB_M.replace("SAVE", save))
        call = f"addpath({_mq(self.dir)}); av_server({_mq(self.dir)})"
        if kind == "octave":
            gui = [] if "octave-cli" in os.path.basename(exe) else ["--no-gui"]
            cmd = [exe, *gui, "--quiet", "--no-window-system", "--eval", call]
        else:
            cmd = [exe, "-batch", call]
        self.logf = open(os.path.join(self.dir, "session.log"), "w")
        group = {"creationflags": subprocess.CREATE_NEW_PROCESS_GROUP} if os.name == "nt" else {"start_new_session": True}
        self.proc = subprocess.Popen(cmd, stdout=self.logf, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, cwd=self.dir, **group)

    def alive(self):
        return self.proc.poll() is None

    def job_dir(self):
        self.n += 1
        d = os.path.join(self.dir, f"job{self.n}")
        os.makedirs(d)
        return d

    def run(self, d, script, timeout, cancelled):
        # Returns (error message or "", printed output); kills the process on timeout or cancel.
        tmp = os.path.join(self.dir, "job.tmp")
        with open(tmp, "w") as fh: fh.write(f"{d}\n{script}")
        os.replace(tmp, os.path.join(self.dir, "job.txt"))
        t0, done = time.monotonic(), os.path.join(d, "done.txt")
        while not os.path.exists(done):
            if not self.alive():
                with open(self.logf.name) as fh: raise RuntimeError(f"{self.kind} session ended: {fh.read()[-800:]}")
            if cancelled() or time.monotonic() - t0 > timeout:
                self.close(kill=True)
                if cancelled(): raise Cancelled()
                raise RuntimeError(f"{self.kind} did not finish within {timeout:g} s (session restarted)")
            time.sleep(0.01)
        with open(done) as fh: err = fh.read()
        log = os.path.join(d, "log.txt")
        out = open(log).read().strip() if os.path.exists(log) else ""
        return err, out

    def kill(self):
        # The matlab command can be a launcher script, so the whole process group/tree is killed, then reaped.
        if os.name == "nt":
            subprocess.run(["taskkill", "/T", "/F", "/PID", str(self.proc.pid)], capture_output=True)
        else:
            import signal
            try: os.killpg(self.proc.pid, signal.SIGKILL)
            except ProcessLookupError: pass
        self.proc.wait()

    def close(self, kill=False):
        if self.alive():
            if kill:
                self.kill()
            else:
                with open(os.path.join(self.dir, "job.txt"), "w") as fh: fh.write("quit")
                try: self.proc.wait(5)
                except subprocess.TimeoutExpired: self.kill()
        self.logf.close()
        shutil.rmtree(self.dir, ignore_errors=True)


_SESSIONS: dict = {}
_SESSION_LOCK = threading.Lock()


def _session(kind, exe):
    s = _SESSIONS.get((kind, exe))
    if s is None or not s.alive():
        s = _SESSIONS[(kind, exe)] = _Session(kind, exe)
    return s


def close_sessions():
    for s in list(_SESSIONS.values()):
        s.close()
    _SESSIONS.clear()


import atexit  # noqa: E402
atexit.register(close_sessions)


class MatlabScript(Modifier):
    title = "MATLAB / Octave script"
    category = "Scripting"
    description = "Run a .m script on every frame using MATLAB (if installed) or free GNU Octave."
    params = [
        Param("file", "file", "", help="A .m script file; leave empty to use the code below"),
        Param("code", "matlab", MATLAB_TEMPLATE),
        Param("engine", "choice", "auto", choices=("auto", "octave", "matlab")),
        Param("timeout", "float", 300.0, min=1.0, label="Timeout (s)"),
        Param("persistent", "bool", True, label="Keep MATLAB running",
              help="Reuse one MATLAB/Octave process for all frames (fast); off = a fresh process per frame"),
    ]

    def label(self):
        return f"MATLAB: {os.path.basename(self.file)}" if self.file else "MATLAB / Octave script"

    def _engine(self):
        engines = find_engines()
        if self.engine == "auto":
            for k in ("matlab", "octave"):
                if k in engines:
                    return k, engines[k]
        elif self.engine in engines:
            return self.engine, engines[self.engine]
        raise RuntimeError("Neither MATLAB nor GNU Octave was found. Install the free Octave with "
                           "'brew install octave', or install MATLAB, then try again.")

    def apply(self, frame, ctx):
        from scipy.io import loadmat, savemat
        kind, exe = self._engine()
        with tempfile.TemporaryDirectory(prefix="atomview_m_") as tmp:
            names = {}  # matlab name -> property name
            variables = {"Position": frame.positions, "Cell": frame.cell, "Origin": frame.origin[None, :],
                         "PBC": np.array(frame.pbc, dtype=float)[None, :], "N": float(frame.n),
                         "Frame": float(ctx.frame_index),
                         "Timestep": float(frame.attributes.get("Timestep", 0) or 0)}
            for name, arr in frame.properties.items():
                a = np.asarray(arr)
                if a.dtype.kind not in "biuf":
                    continue
                mname = _mat_name(name)
                names[mname] = name
                variables[mname] = a.reshape(frame.n, -1).astype(float)
            if self.persistent:
                with _SESSION_LOCK:
                    sess = _session(kind, exe)
                    d = sess.job_dir()
                    savemat(os.path.join(d, "in.mat"), variables)
                    if self.file:
                        script = os.path.abspath(os.path.expanduser(self.file))
                    else:
                        script = os.path.join(d, f"av_s{sess.n}.m")  # unique name: MATLAB caches scripts by name
                        with open(script, "w") as fh: fh.write(self.code)
                    err, self.log = sess.run(d, script, self.timeout, (lambda: ctx.cancelled()) if ctx is not None else (lambda: False))
                    if err: raise RuntimeError(f"{kind}: {err}" + (f"\n{self.log[-600:]}" if self.log else ""))
                    result = loadmat(os.path.join(d, "out.mat"))
                    shutil.rmtree(d, ignore_errors=True)
                return self._merge(frame, result, names, kind)
            inp, outp = os.path.join(tmp, "av_in.mat"), os.path.join(tmp, "av_out.mat")
            savemat(inp, variables)
            if self.file:
                script = os.path.abspath(os.path.expanduser(self.file))
            else:
                script = os.path.join(tmp, "av_user_script.m")
                with open(script, "w") as fh:
                    fh.write(self.code)
            save = f"save('-v7', {_mq(outp)});" if kind == "octave" else f"save({_mq(outp)}, '-v7');"
            wrapper = os.path.join(tmp, "av_wrapper.m")
            with open(wrapper, "w") as fh:
                fh.write(f"load({_mq(inp)});\nrun({_mq(script)});\nclear ans;\n{save}\n")
            if kind == "octave":
                gui_flag = [] if "octave-cli" in os.path.basename(exe) else ["--no-gui"]
                cmd = [exe, *gui_flag, "--quiet", "--no-window-system", "--eval", f"run({_mq(wrapper)})"]
            else:
                cmd = [exe, "-batch", f"run({_mq(wrapper)})"]
            logf = open(os.path.join(tmp, "av_log.txt"), "w+")  # a file, not a pipe: large output cannot block the process
            proc = subprocess.Popen(cmd, stdout=logf, stderr=subprocess.STDOUT, text=True, cwd=tmp, stdin=subprocess.DEVNULL)
            t0 = time.monotonic()
            while proc.poll() is None:  # wait, but stay cancellable and honour the timeout
                if (ctx is not None and ctx.cancelled()) or time.monotonic() - t0 > self.timeout:
                    proc.kill()
                    proc.wait()
                    logf.close()
                    if ctx is not None and ctx.cancelled(): raise Cancelled()
                    raise RuntimeError(f"{kind} did not finish within {self.timeout:g} s")
                time.sleep(0.1)
            logf.seek(0)
            self.log = logf.read().strip()
            logf.close()
            if proc.returncode != 0 or not os.path.exists(outp):
                raise RuntimeError(f"{kind} failed: {self.log[-800:] or 'no output'}")
            result = loadmat(outp)
        return self._merge(frame, result, names, kind)

    def _merge(self, frame, result, names, kind):
        n = frame.n
        skip = {"N", "Frame", "PBC", "Origin"}
        delete = None
        new_props = []
        for key, val in result.items():
            if key.startswith(("__", "av_")) or key in skip or not isinstance(val, np.ndarray) or val.dtype.kind not in "biuf":
                continue
            if key == "Position" and val.shape == (n, 3):
                frame.positions = val.astype(float)
            elif key == "Cell" and val.shape == (3, 3):
                frame.cell = val.astype(float)
            elif key == "DeleteMask" and val.size == n:
                delete = val.reshape(n).astype(bool)
            elif val.ndim == 2 and val.shape[0] == n and n > 1:
                arr = val[:, 0] if val.shape[1] == 1 else val
                prop = names.get(key, key)
                if prop == SELECTION or val.dtype == bool:
                    arr = arr.astype(np.int8) if prop == SELECTION else arr.astype(bool)
                elif np.all(arr == np.round(arr)) and np.asarray(frame.properties.get(prop, arr)).dtype.kind in "iu":
                    arr = arr.astype(np.int64)
                if prop not in frame.properties:
                    new_props.append(prop)
                frame.properties[prop] = arr
            elif val.size == 1 and key != "Timestep":
                frame.attributes[key] = float(val.ravel()[0])
        if delete is not None:
            frame = frame.subset(~delete)
        self.status = f"{kind} ok" + (f"; new: {', '.join(new_props)}" if new_props else "")
        return frame
