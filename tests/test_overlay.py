import numpy as np
import atomview as av
from atomview.expressions import format_overlay
from atomview.render import RenderSettings

def test_format_overlay():
    f=av.Frame(np.zeros((7,3)),attributes={'Timestep':5000,'CNA.counts.FCC':42,'MSD':0.12345,'_hidden':1,'name':'x'})
    assert format_overlay('t = {Timestep*0.002:.1f} ps',f)=='t = 10.0 ps'
    assert format_overlay('FCC {CNA_counts_FCC} of {N}, MSD {MSD:.2e}, frame {Frame}',f,3)=='FCC 42 of 7, MSD 1.23e-01, frame 3'
    assert format_overlay('{max(Timestep, 1)/1000:.0f} k',f)=='5 k'
    assert 'NameError' in format_overlay('{nope}',f) and format_overlay('no fields',f)=='no fields'

def test_overlay_renders_in_animation(tmp_path):
    class S:
        n_frames=3
        def get(s,i): return av.Frame(np.random.default_rng(i).random((30,3))*5,np.eye(3)*5,attributes={'Timestep':i*100,'SourceFrame':i})
    p=av.Pipeline(S()); out=tmp_path/'m.gif'
    av.render_animation(p,str(out),size=(200,150),settings=RenderSettings(overlay='step {Timestep}'))
    import imageio.v2 as iio
    fr=iio.mimread(str(out)); assert len(fr)==3
