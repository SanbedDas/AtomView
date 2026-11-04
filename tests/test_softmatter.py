import numpy as np
from scipy.special import j0
import atomview as av
from atomview import modifiers as m
from atomview.io import mark_2d
from atomview.trajectory import dynamics

class Src:
    def __init__(s,frames): s.fr=frames; s.n_frames=len(frames)
    def get(s,i): return s.fr[i]

def one(f,*mods):
    p=av.Pipeline(Src([f]))
    for x in mods: p.add(x)
    out=p.compute(0); assert not p.errors(), p.errors(); return out

def fcc(n=6,a=1.0):
    b=np.array([[0,0,0],[.5,.5,0],[.5,0,.5],[0,.5,.5]]); g=np.stack(np.meshgrid(*[np.arange(n)]*3,indexing='ij'),-1).reshape(-1,3)
    return ((g[:,None]+b[None]).reshape(-1,3))*a, np.eye(3)*n*a

def tri2d(n=20):
    i,j=np.meshgrid(np.arange(n),np.arange(n),indexing='ij'); x=(i+0.5*(j%2)); y=j*np.sqrt(3)/2
    return mark_2d(av.Frame(np.column_stack([x.ravel(),y.ravel(),np.zeros(x.size)]),np.diag([n,n*np.sqrt(3)/2,1.0])))

def test_structure_factor_gas_and_fcc_peak():
    rng=np.random.default_rng(0); L=20.0
    g=one(av.Frame(rng.uniform(0,L,(8000,3)),np.eye(3)*L),m.StructureFactor(k_max=12,bins=40)).tables['S(k)']
    ok=np.isfinite(g['y']); assert abs(np.nanmean(g['y'][ok][5:])-1)<0.1
    p,c=fcc(); out=one(av.Frame(p,c),m.StructureFactor(k_max=12,bins=200))
    assert abs(out.attributes['StructureFactor.peak_k']-2*np.pi*np.sqrt(3))<0.1  # (111) peak of a unit-lattice FCC

def test_polar_nematic_limits():
    f=tri2d(); n=f.n
    f['theta']=np.zeros(n); out=one(f,m.OrientationalOrder(property='theta'))
    assert np.isclose(out.attributes['Order.polar'],1) and np.isclose(out.attributes['Order.nematic'],1)
    f['theta']=np.where(np.arange(n)%2,0,np.pi); out=one(f,m.OrientationalOrder(property='theta'))
    assert out.attributes['Order.polar']<1e-9 and np.isclose(out.attributes['Order.nematic'],1)
    f['theta']=np.random.default_rng(1).uniform(0,2*np.pi,n); out=one(f,m.OrientationalOrder(property='theta',local_cutoff=1.1))
    assert out.attributes['Order.polar']<0.1 and out.attributes['Order.nematic']<0.15 and out['local_polar_order'].max()<=1+1e-12
    phi=2*np.pi*np.arange(500)/500; v=np.column_stack([np.cos(phi),np.sin(phi),np.zeros(500)])  # planar disk in 3D: S = 1/4
    out=one(av.Frame(np.random.default_rng(3).random((500,3))*10,np.eye(3)*10,properties={'u':v}),m.OrientationalOrder(property='u'))
    assert np.isclose(out.attributes['Order.nematic'],0.25,atol=1e-12) and out.attributes['Order.polar']<1e-12

def test_g6_perfect_lattice():
    t=one(tri2d(),m.BondOrderCorrelation(r_max=5,bins=50)).tables['g6(r)']
    assert np.allclose(t['y'][np.isfinite(t['y'])],1,atol=1e-6)

def ballistic(dim,F=30,N=400,v=0.07,L=50.0):
    rng=np.random.default_rng(4); u=rng.normal(size=(N,dim)); u/=np.linalg.norm(u,axis=1)[:,None]  # same speed, random direction
    p0=rng.uniform(0,L,(N,dim)); frames=[]
    for t in range(F):
        p=(p0+u*v*t)%L; p=np.column_stack([p,np.zeros(N)]) if dim==2 else p
        f=av.Frame(p,np.diag([L,L,L if dim==3 else 1.0]),pbc=(True,True,dim==3),properties={'id':np.arange(1,N+1),'velocity':np.column_stack([u*v,np.zeros(N)]) if dim==2 else u*v},attributes={'Timestep':t,**({'Dimensions':2} if dim==2 else {})})
        frames.append(f)
    return av.Pipeline(Src(frames)),v

def test_dynamics_ballistic_exact():
    for dim in (3,2):
        p,v=ballistic(dim); k=2.0; r=dynamics(p,k=k)
        t=r['MSD(t)']['x']; kv=k*v*t
        assert np.allclose(r['MSD(t)']['y'],(v*t)**2,rtol=1e-9)
        assert np.allclose(r['alpha2(t)']['y'][1:],dim/(dim+2)-1,atol=1e-9)  # delta distribution of |dr|
        assert np.allclose(r[f'Fs(k={k:g},t)']['y'],j0(kv) if dim==2 else np.sinc(kv/np.pi),atol=1e-9)
        assert np.allclose(r['VACF(t)']['y'],1,atol=1e-9) and r['VACF(t)']['source']=='velocity property'

def test_vacf_uncorrelated_steps():
    rng=np.random.default_rng(5); N,F,L=3000,40,100.0; P=np.cumsum(rng.normal(0,.2,(F,N,3)),axis=0)+50
    fr=[av.Frame(P[t]%L,np.eye(3)*L,properties={'id':np.arange(1,N+1)}) for t in range(F)]
    c=dynamics(av.Pipeline(Src(fr)))['VACF(t)']
    assert c['source']=='finite differences' and abs(c['y'][1])<0.03 and abs(c['y'][5])<0.03

def test_trails(tmp_path):
    p,v=ballistic(2,F=8)
    p.add(m.Trails(length=5)); out=p.compute(7); tr=out.attributes['_trails']
    assert tr['path'].shape==(6,400,3) and np.allclose(np.diff(tr['path'],axis=0)[...,:2],(p.compute(0)['velocity'][:,:2])[None],atol=1e-9)
    av.render_image(out,str(tmp_path/'t.png'),size=(300,200)); assert (tmp_path/'t.png').stat().st_size>1000
