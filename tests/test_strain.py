import numpy as np
import atomview as av
from atomview import modifiers as m
from atomview.io import mark_2d

def fcc(n=5,a=3.6):
    b=np.array([[0,0,0],[.5,.5,0],[.5,0,.5],[0,.5,.5]]); g=np.stack(np.meshgrid(*[np.arange(n)]*3,indexing='ij'),-1).reshape(-1,3)
    return ((g[:,None]+b[None]).reshape(-1,3))*a, np.eye(3)*n*a

class Two:
    def __init__(s,f0,f1): s.f=[f0,f1]; s.n_frames=2
    def get(s,i): return s.f[i]

def deformed(p,c,F,order=None):
    o=np.arange(len(p)) if order is None else order
    return av.Frame(p[o]@F.T,c@F.T,properties={'id':np.arange(1,len(p)+1)[o]})

def run(f0,f1,cut=3.0,**kw):
    pp=av.Pipeline(Two(f0,f1)); pp.add(m.AtomicStrain(cutoff=cut,output_tensor=True,output_gradient=True,**kw)); out=pp.compute(1)
    assert not pp.errors(), pp.errors(); return out

def test_simple_shear_exact():
    p,c=fcc(); g=0.05; F=np.eye(3); F[0,1]=g
    out=run(av.Frame(p,c,properties={'id':np.arange(1,len(p)+1)}),deformed(p,c,F))
    assert np.allclose(out['shear_strain'],np.sqrt(g**2/4+g**4/12),atol=1e-12)
    assert np.allclose(out['volumetric_strain'],g**2/6,atol=1e-12) and np.allclose(out['D2min'],0,atol=1e-20)
    assert np.allclose(out['deformation_gradient'],F.ravel(),atol=1e-12) and out.selection.sum()==0
    assert np.allclose(out['strain_tensor'],[0,g**2/2,0,g/2,0,0],atol=1e-12)

def test_expansion_and_shuffled_order():
    p,c=fcc(); F=np.eye(3)*1.01; order=np.random.default_rng(0).permutation(len(p))
    out=run(av.Frame(p,c,properties={'id':np.arange(1,len(p)+1)}),deformed(p,c,F,order))
    assert np.allclose(out['shear_strain'],0,atol=1e-12) and np.allclose(out['volumetric_strain'],(1.01**2-1)/2,atol=1e-12)

def test_nonaffine_move_is_local():
    p,c=fcc(); q=p.copy(); q[0]+=[0.2,-0.1,0.15]
    out=run(av.Frame(p,c,properties={'id':np.arange(1,len(p)+1)}),av.Frame(q,c,properties={'id':np.arange(1,len(p)+1)}))
    near=np.linalg.norm(av.Frame(p,c).minimum_image(p-p[0]),axis=1)<3.0
    assert out['D2min'][0]>0.05 and (out['D2min'][near]>0).all() and np.allclose(out['D2min'][~near],0,atol=1e-20)

def test_2d_shear():
    i,j=np.meshgrid(np.arange(16),np.arange(16),indexing='ij'); x=(i+0.5*(j%2)); y=j*np.sqrt(3)/2
    p=np.column_stack([x.ravel(),y.ravel(),np.zeros(x.size)]); c=np.diag([16,16*np.sqrt(3)/2,1.0]); g=0.04
    F=np.eye(3); F[0,1]=g
    f0=mark_2d(av.Frame(p,c,properties={'id':np.arange(1,len(p)+1)})); f1=mark_2d(deformed(p,c,F))
    out=run(f0,f1,cut=1.2)
    E=np.array([[0,g/2],[g/2,g**2/2]]); assert np.allclose(out['shear_strain'],np.sqrt(E[0,1]**2+(E[0,0]-E[1,1])**2/2),atol=1e-12)
