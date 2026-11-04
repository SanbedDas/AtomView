import numpy as np
import pytest
from scipy.io import savemat
import atomview as av
from atomview import modifiers as m
from atomview.io import mark_2d, read_mat

def tri(nx,ny,a=1.0):
    # periodic triangular lattice, rows offset by a/2, row spacing a*sqrt(3)/2
    i,j=np.meshgrid(np.arange(nx),np.arange(ny),indexing='ij'); x=(i+0.5*(j%2))*a; y=j*a*np.sqrt(3)/2
    return np.column_stack([x.ravel(),y.ravel()]), nx*a, ny*a*np.sqrt(3)/2

def frame2d(p2,Lx,Ly):
    return mark_2d(av.Frame(np.column_stack([p2,np.zeros(len(p2))]),np.diag([Lx,Ly,1.0])))

class Src:
    def __init__(s,f): s.f=f; s.n_frames=1
    def get(s,i): return s.f

def test_triangular_lattice_order_and_coordination():
    f=frame2d(*tri(20,20))
    assert f.attributes['Dimensions']==2 and f.pbc==(True,True,False)
    p=av.Pipeline(Src(f)); p.add(m.Hexatic()); p.add(m.Coordination(cutoff=1.2))
    out=p.compute(0); assert not p.errors(), p.errors()
    assert np.allclose(out['psi6_abs'],1,atol=1e-4) and out.attributes['Hexatic.global_psi6']>0.999
    assert (out['coordination']==6).all()

def test_ideal_gas_rdf_is_one():
    rng=np.random.default_rng(1); L=60.0
    f=frame2d(rng.uniform(0,L,(20000,2)),L,L)
    p=av.Pipeline(Src(f)); p.add(m.Coordination(cutoff=5.0,bins=50))
    g=p.compute(0).tables['Radial distribution function']['y']
    assert abs(g[10:].mean()-1)<0.03

def test_3d_only_modifiers_refuse_2d():
    f=frame2d(*tri(10,10)); p=av.Pipeline(Src(f)); p.add(m.CommonNeighborAnalysis()); p.add(m.Steinhardt())
    p.compute(0); assert [e for _,e in p.errors() if 'Hexatic' in e]==[e for _,e in p.errors()] and len(p.errors())==2

def test_mat_2d_trajectory(tmp_path):
    p2,Lx,Ly=tri(8,8); traj=np.stack([p2+0.01*k for k in range(5)],axis=2)  # N x 2 x F
    savemat(tmp_path/'c.mat',{'pos':traj,'BoxSize':[Lx,Ly],'theta':np.zeros((len(p2),5))})
    fr=read_mat(str(tmp_path/'c.mat'))
    assert len(fr)==5 and fr[0].attributes['Dimensions']==2 and fr[0].pbc==(True,True,False)
    assert np.allclose(np.diag(fr[0].cell),[Lx,Ly,1]) and fr[4].positions[0,0]==pytest.approx(0.04)
    assert fr[2]['theta'].shape==(len(p2),)
    t=av.load(str(tmp_path/'c.mat')); assert t.n_frames==5 and t.get(1).attributes['Dimensions']==2

def test_flat_3d_input_detected_as_2d(tmp_path):
    p2,Lx,Ly=tri(6,6); fn=tmp_path/'flat.xyz'
    with open(fn,'w') as fh:
        fh.write(f'{len(p2)}\nLattice="{Lx} 0 0 0 {Ly} 0 0 0 1" Properties=species:S:1:pos:R:3\n')
        for x,y in p2: fh.write(f'A {x} {y} 0.0\n')
    f=av.load(str(fn)).get(0); assert f.attributes['Dimensions']==2 and not f.pbc[2]

def test_render_2d(tmp_path):
    f=frame2d(*tri(12,12)); out=tmp_path/'t.png'
    av.render_image(f,str(out),size=(400,300)); assert out.stat().st_size>1000

def test_vectors_from_angle_and_render(tmp_path):
    from atomview.modifiers import vector_field
    p2,Lx,Ly=tri(6,6); f=frame2d(p2,Lx,Ly); f['theta']=np.full(f.n,np.pi/2); f['v2']=np.tile([3.0,4.0],(f.n,1))
    assert np.allclose(vector_field(f,'theta'),[0,1,0]) and np.allclose(vector_field(f,'v2'),[3,4,0])
    p=av.Pipeline(Src(f)); p.add(m.Vectors(property='theta')); p.add(m.DeleteSelected()); out=p.compute(0)
    assert not p.errors() and out.attributes['_vectors'][0]['property']=='theta'
    av.render_image(out,str(tmp_path/'v.png'),size=(300,200)); assert (tmp_path/'v.png').stat().st_size>1000

def test_set_cell_centered_box():
    p,Lx,Ly=tri(10,10); f=frame2d(p-5.0,Lx,Ly)
    g=m.SetCell(lx=20,ly=40,pbc="y").apply(f.copy(),None)
    assert np.allclose(np.diag(g.cell)[:2],[20,40]) and np.allclose(g.origin[:2],[-10,-20])
    assert g.pbc==(False,True,False) and g.cell[2,2]==f.cell[2,2] and np.array_equal(g.positions,f.positions)
