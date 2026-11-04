import numpy as np
import pytest
import atomview as av
from atomview.trajectory import msd

class Traj:
    def __init__(s,P,L,dim=3,ts=None):
        s.P,s.L,s.dim,s.ts=P,L,dim,ts; s.n_frames=len(P)
    def get(s,k):
        p=s.P[k]%s.L  # stored wrapped, as simulations write them
        a={'Timestep':int(s.ts[k])} if s.ts is not None else {}
        if s.dim==2: a['Dimensions']=2
        return av.Frame(p,np.diag([s.L,s.L,s.L if s.dim==3 else 1.0]),pbc=(True,True,s.dim==3),
                        properties={'id':np.arange(1,len(p)+1)},attributes=a)

def test_ballistic_exact_through_periodic_wrapping():
    rng=np.random.default_rng(0); N,F,L=200,40,10.0; v=rng.normal(0,0.3,(N,3))
    P=rng.uniform(0,L,(N,3))[None]+np.arange(F)[:,None,None]*v[None]  # crosses the box many times
    r=msd(av.Pipeline(Traj(P,L)),dt=0.5)
    t=np.arange(F)*0.5
    assert np.allclose(r['y'],np.mean((v**2).sum(1))*(t/0.5)**2,rtol=1e-9) and np.allclose(r['x'],t)

def test_2d_random_walk_diffusion():
    rng=np.random.default_rng(1); N,F,L,s=2000,200,50.0,0.1
    P=np.cumsum(np.concatenate([rng.uniform(0,L,(1,N,2)),rng.normal(0,s,(F-1,N,2))]),axis=0)
    P=np.concatenate([P,np.zeros((F,N,1))],axis=2)
    r=msd(av.Pipeline(Traj(P,L,dim=2,ts=np.arange(F)*100)),fit_from=0.0)
    assert r['dim']==2 and r['x'][1]==100
    D_true=s**2/2/100  # per unit time, timestep 100 per frame
    assert abs(r['D']/D_true-1)<0.05, (r['D'],D_true)

def test_strided_frames_and_errors():
    P=np.cumsum(np.random.default_rng(2).normal(0,.1,(30,50,3)),axis=0)+5
    r=msd(av.Pipeline(Traj(P,10.0)),frames=range(0,30,3),dt=2.0)
    assert np.allclose(r['x'][:3],[0,6,12])
    with pytest.raises(ValueError): msd(av.Pipeline(Traj(P,10.0)),frames=[0,1,3])

def test_fft_msd_matches_direct_sum():
    from atomview.trajectory import windowed_msd
    r=np.cumsum(np.random.default_rng(4).normal(0,1,(25,7,3)),axis=0)
    direct=[np.mean([((r[t+m]-r[t])**2).sum(-1).mean() for t in range(25-m)]) for m in range(25)]
    assert np.allclose(windowed_msd(r),direct,rtol=1e-10,atol=1e-12)

def test_particles_missing_in_some_frames_are_dropped():
    class Leaky(Traj):
        def get(s,k):
            f=super().get(k); return f.subset(f['id']!=(k%5)+1) if k else f  # one particle missing per frame
    P=np.cumsum(np.random.default_rng(5).normal(0,.05,(12,40,3)),axis=0)+5
    full=msd(av.Pipeline(Traj(P[:,5:],10.0)))
    r=msd(av.Pipeline(Leaky(P,10.0)))
    assert r['dropped']==5 and r['particles']==35 and np.allclose(r['y'],full['y'])
