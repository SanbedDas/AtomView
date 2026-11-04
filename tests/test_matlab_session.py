import os, time, threading
import numpy as np
import pytest
import atomview as av
from atomview import modifiers as m
from atomview.pipeline import Cancelled

pytestmark = pytest.mark.skipif(not m.find_engines(), reason="MATLAB/Octave not installed")
DATA = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples", "data", "cu_void.dump")

def pipe(**kw):
    p = av.Pipeline(av.load(DATA)); mod = p.add(m.MatlabScript(**kw)); return p, mod

def test_persistent_equals_one_shot_and_is_isolated():
    code = "if exist('leak','var'), leaked = 1; else, leaked = 0; end\nleak = 1;\nr = sqrt(sum(Position.^2, 2));\nselection = r < 20;\nMeanR = mean(r);\n"
    a, _ = pipe(code=code, persistent=True); b, _ = pipe(code=code, persistent=False)
    for f in (0, 1, 2):
        x, y = a.compute(f), b.compute(f)
        assert not a.errors() and not b.errors()
        assert np.array_equal(x["r"], y["r"]) and np.array_equal(x.selection, y.selection)
        assert x.attributes["MeanR"] == y.attributes["MeanR"] and x.attributes["leaked"] == 0

def test_error_keeps_session_and_file_edits_are_seen(tmp_path):
    fn = tmp_path / "s.m"; fn.write_text("v = 1;\n")
    p, mod = pipe(file=str(fn))
    assert p.compute(0).attributes["v"] == 1
    fn.write_text("v = 2;\n"); os.utime(fn, (time.time() + 2, time.time() + 2))
    p.invalidate(); assert p.compute(0).attributes["v"] == 2
    fn.write_text("v = undefined_function_xyz(3);\n"); os.utime(fn, (time.time() + 4, time.time() + 4))
    p.invalidate(); p.compute(0); assert "undefined_function_xyz" in mod.error
    fn.write_text("v = 3;\n"); os.utime(fn, (time.time() + 6, time.time() + 6))
    p.invalidate(); assert p.compute(0).attributes["v"] == 3

def test_cancel_kills_and_next_run_restarts():
    p, mod = pipe(code="pause(60);\n"); p.strict = True
    threading.Timer(2.0, p.cancel.set).start(); t = time.perf_counter()
    with pytest.raises(Cancelled): p.compute(0)
    assert time.perf_counter() - t < 5
    p.cancel.clear(); p.strict = False; mod.values["code"] = "ok = 1;\n"; p.invalidate()
    assert p.compute(0).attributes["ok"] == 1
