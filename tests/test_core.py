import os
import subprocess
import sys

import numpy as np
import pytest

import atomview as av
from atomview import modifiers as m
from atomview.io import Trajectory, write

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "examples", "data")


@pytest.fixture(scope="session", autouse=True)
def samples():
    if not os.path.exists(os.path.join(DATA, "cu_void.dump")):
        subprocess.run([sys.executable, os.path.join(ROOT, "examples", "make_samples.py")], check=True)


def test_dump_reader():
    t = Trajectory(os.path.join(DATA, "cu_void.dump"))
    assert t.n_frames == 20
    f = t.get(3)
    assert f.attributes["Timestep"] == 3000
    assert "c_pe" in f.properties and f.properties["type"].dtype.kind == "i"
    assert np.allclose(np.diag(f.cell), 36.15)


def test_xyz_reader_elements():
    f = av.load(os.path.join(DATA, "fe_bcc.xyz")).get(0)
    assert f.type_names == {1: "Fe"} and f.n == 1024


def test_cna_fcc_and_bcc():
    p = av.Pipeline(av.load(os.path.join(DATA, "cu_void.dump")))
    p.add(m.CommonNeighborAnalysis())
    out = p.compute(0)
    st = out.properties["structure_type"]
    assert (st == 1).mean() > 0.8  # mostly FCC; surface of the void is "Other"
    q = av.Pipeline(av.load(os.path.join(DATA, "fe_bcc.xyz")))
    q.add(m.CommonNeighborAnalysis())
    assert (q.compute(0).properties["structure_type"] == 3).mean() > 0.95


def test_hcp_detected():
    from ase.build import bulk
    from atomview.io import from_ase
    f = from_ase(bulk("Mg", "hcp", a=3.2, c=5.2).repeat((5, 5, 4)))
    st = m.adaptive_cna(f)
    assert (st == 2).all()


def test_expression_select_and_delete():
    p = av.Pipeline(av.load(os.path.join(DATA, "cu_void.dump")))
    p.add(m.ExpressionSelect(expression="type == 2 && c_pe > -3.5"))
    sel = p.compute(0).selection.sum()
    assert sel > 0
    p.add(m.DeleteSelected())
    assert p.compute(0).n == p.source.get(0).n - sel


def test_python_script_and_function():
    p = av.Pipeline(av.load(os.path.join(DATA, "cu_void.dump")))
    p.add(m.PythonScript(code="def modify(frame, data):\n    data['double_pe'] = 2 * data['c_pe']\n    print('hi')\n"))
    p.add(lambda frame, data: data.attributes.update(FrameSeen=frame))
    out = p.compute(2)
    assert np.allclose(out["double_pe"], 2 * out["c_pe"])
    assert out.attributes["FrameSeen"] == 2
    assert p.modifiers[0].log.strip() == "hi"


def test_analysis_modifiers_run():
    p = av.Pipeline(av.load(os.path.join(DATA, "cu_void.dump")))
    for mod in (m.Coordination(cutoff=3.0), m.CentroSymmetry(), m.ClusterAnalysis(cutoff=3.0),
                m.CreateBonds(cutoff=2.8), m.Steinhardt(), m.VoronoiAnalysis(), m.Displacement(),
                m.ColorCoding(property="c_pe"), m.Histogram(property="c_pe"), m.Slice(width=10),
                m.Replicate(nx=2, ny=1, nz=1), m.WrapPeriodic()):
        p.add(mod)
    out = p.compute(5)
    assert not p.errors(), p.errors()
    assert "Radial distribution function" in out.tables
    assert 10 < out.properties["coordination"].mean() <= 12
    assert out.bonds is not None and len(out.bonds)


def test_session_roundtrip(tmp_path):
    p = av.Pipeline(av.load(os.path.join(DATA, "cu_void.dump")))
    p.add(m.ColorCoding(property="c_pe", colormap="plasma"))
    path = tmp_path / "s.json"
    p.save(path)
    q = av.Pipeline.load(path)
    assert q.modifiers[0].colormap == "plasma"


def test_writers_roundtrip(tmp_path):
    f = av.load(os.path.join(DATA, "cu_void.dump")).get(0)
    for name in ("o.xyz", "o.dump"):
        write(f, tmp_path / name)
        g = av.load(str(tmp_path / name)).get(0)
        assert g.n == f.n and np.allclose(g.positions, f.positions, atol=1e-4)


@pytest.mark.skipif(not m.find_engines(), reason="MATLAB/Octave not installed")
def test_matlab_script():
    p = av.Pipeline(av.load(os.path.join(DATA, "cu_void.dump")))
    p.add(m.MatlabScript(code="twice = 2 * c_pe;\nselection = Position(:,3) > 10;\nAvgPE = mean(c_pe);\n"))
    out = p.compute(0)
    assert not p.errors(), p.errors()
    assert np.allclose(out["twice"], 2 * out["c_pe"])
    assert abs(out.attributes["AvgPE"] - out["c_pe"].mean()) < 1e-9
