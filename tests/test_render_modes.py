import numpy as np
import atomview as av
from atomview.render import RenderSettings

def cover(fn):
    import imageio.v2 as iio
    a=iio.imread(fn)[...,:3].astype(int); return (np.abs(a-255).sum(-1)>30).mean()

def test_imposter_matches_mesh_coverage(tmp_path):
    rng=np.random.default_rng(0); p=rng.uniform(0,40,(3000,3))
    f=av.Frame(p,np.eye(3)*40,properties={'type':np.ones(3000,int),'radius':np.full(3000,0.8)})
    cam=[(110,-90,95),(20,20,20),(0,0,1)]; out={}
    for mode in ('mesh','imposter'):
        fn=str(tmp_path/f'{mode}.png'); s=RenderSettings(render_mode=mode,show_cell=False,show_legend=False,show_axes=False)
        av.render_image(f,fn,size=(600,450),settings=s,camera=cam); out[mode]=cover(fn)
    assert out['mesh']>0.05 and abs(out['imposter']-out['mesh'])/out['mesh']<0.1, out

def test_auto_threshold():
    from atomview.render import use_imposters
    s=RenderSettings(); small=av.Frame(np.zeros((10,3))); big=av.Frame(np.zeros((40000,3)))
    assert not use_imposters(small,s) and use_imposters(big,s) and not use_imposters(big,RenderSettings(render_mode='mesh'))
