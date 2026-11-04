import threading, time
import numpy as np
import pytest
import atomview as av
from atomview import modifiers as m
from atomview.pipeline import Cancelled

class S:
    n_frames=5
    def get(s,i): return av.Frame(np.random.default_rng(i).random((50,3))*5,np.eye(3)*5,properties={'id':np.arange(1,51)})

def slow(frame,data): time.sleep(0.3); return data

def test_cancel_between_modifiers_and_recover():
    p=av.Pipeline(S())
    for _ in range(10): p.add(slow)
    msgs=[]; p.progress=msgs.append
    threading.Timer(0.5,p.cancel.set).start(); t=time.perf_counter()
    with pytest.raises(Cancelled): p.compute(0)
    assert time.perf_counter()-t<1.2 and not p._cache and 1<=len(msgs)<=3
    p.cancel.clear(); p.modifiers=p.modifiers[:1]; p.invalidate(); assert p.compute(0).n==50

def test_msd_progress_cancel():
    from atomview.trajectory import msd
    with pytest.raises(Cancelled): msd(av.Pipeline(S()),progress=lambda k,n: k<2)

@pytest.mark.skipif(not m.find_engines(),reason="MATLAB/Octave not installed")
def test_matlab_process_is_killed_on_cancel():
    p=av.Pipeline(S()); p.add(m.MatlabScript(code="pause(60);\n"))
    p.strict=True; threading.Timer(3.0,p.cancel.set).start(); t=time.perf_counter()
    with pytest.raises(Cancelled): p.compute(0)
    assert time.perf_counter()-t<6

def test_script_output_capture_is_per_thread(capsys):
    import io
    from atomview.modifiers.scripting import capture_stdout
    started=threading.Event(); buf=io.StringIO()
    def worker():
        with capture_stdout(buf):
            print("from script"); started.set(); time.sleep(0.3)
    th=threading.Thread(target=worker); th.start(); started.wait()
    print("from main"); th.join()
    assert buf.getvalue()=="from script\n" and "from main" in capsys.readouterr().out
