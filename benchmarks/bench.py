# Times load, CNA and offscreen render vs particle count; writes a timestamped log.
import os, sys, time, json, resource
import numpy as np
import atomview as av
from atomview import modifiers as m
from atomview.io import write_lammps_dump
from atomview.render import RenderSettings

def fcc(n,a=3.615):
    b=np.array([[0,0,0],[.5,.5,0],[.5,0,.5],[0,.5,.5]]); g=np.stack(np.meshgrid(*[np.arange(n)]*3,indexing='ij'),-1).reshape(-1,3)
    return ((g[:,None]+b[None]).reshape(-1,3))*a, n*a

def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(2**30 if sys.platform=='darwin' else 2**20)

def run(sizes,mode,out):
    res=[]
    for n in sizes:
        p,L=fcc(n); p=p+np.random.default_rng(0).normal(0,.05,p.shape)
        f=av.Frame(p,np.eye(3)*L,properties={'id':np.arange(1,len(p)+1),'type':np.ones(len(p),int)},attributes={'Timestep':0})
        fn=os.path.join(out,f'fcc{n}.dump'); write_lammps_dump(f,fn)
        t=time.perf_counter(); tr=av.load(fn); fr=tr.get(0); tl=time.perf_counter()-t
        t=time.perf_counter(); m.adaptive_cna(fr); tc=time.perf_counter()-t
        s=RenderSettings(render_mode=mode) if 'render_mode' in RenderSettings.__dataclass_fields__ else RenderSettings()
        t=time.perf_counter(); av.render_image(fr,os.path.join(out,f'fcc{n}_{mode}.png'),size=(1600,1200),settings=s); tr_=time.perf_counter()-t
        r=dict(N=len(p),load_s=round(tl,3),cna_s=round(tc,3),render_s=round(tr_,3),mode=mode,peak_rss_gb=round(rss(),2)); res.append(r); print(r,flush=True)
    return res

if __name__=='__main__':
    out=sys.argv[1] if len(sys.argv)>1 else 'build/bench'; mode=sys.argv[2] if len(sys.argv)>2 else 'auto'
    sizes=[int(x) for x in (sys.argv[3] if len(sys.argv)>3 else '14,29,63').split(',')]  # 4n^3 atoms: 11k, 98k, 1.0M
    os.makedirs(out,exist_ok=True); res=run(sizes,mode,out)
    with open(os.path.join(out,time.strftime('bench_%Y%m%d_%H%M%S.json')),'w') as fh: json.dump(dict(platform=sys.platform,python=sys.version.split()[0],sizes=sizes,mode=mode,results=res),fh,indent=1)
