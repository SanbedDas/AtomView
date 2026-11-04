import os, subprocess, sys, textwrap
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

PLUG='''
import numpy as np
from atomview.pipeline import Modifier, Param
class {cls}(Modifier):
    title="{title}"; category="Plugins"; params=[Param("factor","float",2.0)]
    def apply(self,frame,ctx): frame["scaled_x"]=self.factor*frame.positions[:,0]; return frame
'''

def run(code,env):
    return subprocess.run([sys.executable,"-c",textwrap.dedent(code)],capture_output=True,text=True,env={**os.environ,**env},cwd=ROOT,timeout=120)

def test_folder_and_entry_point_plugins(tmp_path):
    d=tmp_path/"plugdir"; d.mkdir(); (d/"scale.py").write_text(PLUG.format(cls="ScaleX",title="Scale x"))
    (d/"broken.py").write_text("raise RuntimeError('boom')")
    site=tmp_path/"site"; site.mkdir(); (site/"epplug.py").write_text(PLUG.format(cls="EPScale",title="EP scale"))
    di=site/"epplug-0.1.dist-info"; di.mkdir(); (di/"METADATA").write_text("Metadata-Version: 2.1\nName: epplug\nVersion: 0.1\n")
    (di/"entry_points.txt").write_text("[atomview.modifiers]\nepscale = epplug\n")
    r=run('''
        import numpy as np, atomview as av
        from atomview.plugins import load_plugins, ERRORS
        from atomview.pipeline import REGISTRY, Pipeline
        load_plugins()
        assert "ScaleX" in REGISTRY and "EPScale" in REGISTRY, list(REGISTRY)
        assert len(ERRORS)==1 and "boom" in ERRORS[0][1]
        f=av.Frame(np.ones((5,3)))
        class S:
            n_frames=1
            def get(s,i): return f
        p=Pipeline(S()); p.add(REGISTRY["ScaleX"](factor=3.0)); assert np.allclose(p.compute(0)["scaled_x"],3)
        print("ok")
    ''',{"ATOMVIEW_PLUGIN_PATH":str(d),"PYTHONPATH":str(site)})
    assert r.stdout.strip()=="ok", r.stdout+r.stderr
    r=subprocess.run([sys.executable,"-m","atomview","plugins"],capture_output=True,text=True,env={**os.environ,"ATOMVIEW_PLUGIN_PATH":str(d),"PYTHONPATH":str(site)},cwd=ROOT,timeout=120)
    assert "Scale x" in r.stdout and "EP scale" in r.stdout and "boom" in r.stderr and r.returncode==1
