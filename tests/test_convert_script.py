import os, subprocess, sys
import numpy as np
import pytest
import atomview as av
from atomview.io import read_mat, write_mat, _savemat73
from atomview.pipeline import parse_modifier
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DUMP=os.path.join(ROOT,"examples","data","cu_void.dump")

def cli(*a): return subprocess.run([sys.executable,"-m","atomview",*a],capture_output=True,text=True,cwd=ROOT,timeout=300)

def test_convert_roundtrip_with_modifier(tmp_path):
    out=str(tmp_path/"c.mat"); r=cli("convert",DUMP,out,"-m","PolyhedralTemplateMatching()","-m","ComputeProperty(output='h', expression='z*2')","--frames","0:20:5")
    assert r.returncode==0, r.stderr
    p=av.Pipeline(av.load(DUMP)); p.add(parse_modifier("PolyhedralTemplateMatching()")); p.add(parse_modifier("ComputeProperty(output='h', expression='z*2')"))
    back=read_mat(out); assert len(back)==4
    for k,f in zip((0,5,10,15),back):
        ref=p.compute(k)
        assert np.array_equal(f.positions,ref.positions) and np.array_equal(f['structure_type'],ref['structure_type']) and np.allclose(f['h'],ref['h'])
        assert f.attributes['PTM_counts_FCC']==ref.attributes['PTM.counts.FCC'] and f.attributes['Timestep']==ref.attributes['Timestep']

def test_parse_modifier_is_literal_only():
    m=parse_modifier("ColorCoding(property='c_pe', colormap='plasma', start=-1.5)"); assert m.colormap=='plasma' and m.start==-1.5
    assert type(parse_modifier("Polyhedral template matching")).__name__=="PolyhedralTemplateMatching"
    for bad in ("ColorCoding(property=__import__('os').system('echo hi'))","os.system('x')","ColorCoding('c_pe')","Nope()"):
        with pytest.raises((ValueError,TypeError)): parse_modifier(bad)

def test_script_command(tmp_path):
    s=tmp_path/"s.py"; s.write_text("import sys, atomview\nf=atomview.load(sys.argv[1]).get(0)\nprint('N', f.n, sys.argv[2])\n")
    r=cli("script",str(s),DUMP,"--flag"); assert r.returncode==0 and r.stdout.strip()=="N 3921 --flag", r.stdout+r.stderr

def test_v73_writer_reads_back(tmp_path):
    f=av.load(DUMP).get(3); fn=str(tmp_path/"big.mat")
    _savemat73(fn,{"Position":f.positions,"Cell":f.cell,"type":f['type'],"c_pe":f['c_pe']})
    g=read_mat(fn)[0]; assert np.array_equal(g.positions,f.positions) and np.array_equal(g['c_pe'],f['c_pe'])
    with open(fn,"rb") as fh: assert fh.read(19)==b"MATLAB 7.3 MAT-file"
