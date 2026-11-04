import numpy as np
import pytest
from ase.build import bulk
from scipy.spatial.transform import Rotation
import atomview as av
from atomview import modifiers as m
from atomview.io import from_ase
from atomview.ptm import ptm

ALL=("fcc","hcp","bcc","ico","sc","dcub","dhex")

@pytest.mark.parametrize("name,kind,a,rep,want",[("Cu","fcc",3.6,5,1),("Mg","hcp",3.2,(5,5,4),2),("Fe","bcc",2.87,6,3),("Po","sc",3.35,6,5),("Si","diamond",5.43,4,6)])
def test_perfect_lattices(name,kind,a,rep,want):
    at=bulk(name,kind,a=a,c=5.2,cubic=kind!="hcp") if kind=="hcp" else bulk(name,kind,a=a,cubic=True)
    r=ptm(from_ase(at.repeat(rep)),ALL)
    assert (r['type']==want).all() and np.nanmax(r['rmsd'])<3e-3
    if kind!="hcp": assert np.allclose(r['lattice_constant'],a,rtol=1e-6)

def test_orientation_recovered_up_to_cubic_symmetry():
    f=from_ase(bulk("Cu","fcc",a=3.6,cubic=True).repeat(6)); R=Rotation.from_euler('zyx',[23,41,-17],degrees=True).as_matrix()
    g=av.Frame(f.positions@R.T,f.cell@R.T,pbc=f.pbc)  # rigidly rotated crystal and cell
    q=ptm(g)['q']; Rp=Rotation.from_quat(np.roll(q,-1,axis=1)).as_matrix()  # PTM (w,x,y,z) -> scipy (x,y,z,w)
    M=np.einsum('ij,nik->njk',R,Rp)  # R^T Rp must be a cubic point-group operation: entries 0 or +-1
    assert np.allclose(np.abs(M).round(),np.abs(M),atol=1e-6) and np.allclose(np.abs(M).sum(1),1,atol=1e-6)

def test_thermal_noise_and_rmsd_cutoff_vs_cna():
    f=av.load('examples/data/cu_void.dump').get(10)
    p=av.Pipeline(type("S",(),{"n_frames":1,"get":lambda s,i: f})())
    p.add(m.PolyhedralTemplateMatching(rmsd_cutoff=0)); out=p.compute(0); assert not p.errors()
    ptm_fcc=(out['structure_type']==1); cna_fcc=(m.adaptive_cna(f)==1)
    assert ptm_fcc.mean()>=cna_fcc.mean() and (ptm_fcc|~cna_fcc).mean()>0.99  # PTM finds every CNA-FCC atom (and more)
    assert 'PTM RMSD' in out.tables and out['orientation'].shape==(f.n,4)
    p.modifiers[0].values['rmsd_cutoff']=0.05; p.invalidate(); st=p.compute(0)['structure_type']
    assert (st==0).sum()>(out['structure_type']==0).sum()

def test_alloy_ordering_L12():
    at=bulk("Cu","fcc",a=3.6,cubic=True).repeat(4); f=from_ase(at)
    fc=f.fractional()*4 % 1; corner=np.all(np.isclose(fc,0,atol=1e-6)|np.isclose(fc,1,atol=1e-6),axis=1)
    types=np.where(corner,2,1)  # Cu3Au: Au on the cube corners
    r=ptm(f,("fcc",),numbers=types)
    assert set(np.unique(r['alloy'][types==1]))=={3} and set(np.unique(r['alloy'][types==2]))=={4}  # L12 Cu / L12 Au
