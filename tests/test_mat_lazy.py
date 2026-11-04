import numpy as np
import h5py
import atomview as av
from atomview.io import MatReader, read_mat, _H5Lazy

def test_v73_frames_read_on_demand(tmp_path):
    rng=np.random.default_rng(3); N,F=70000,6; P=rng.random((N,3,F)); V=rng.random((N,F))
    fn=str(tmp_path/'t.mat')
    with h5py.File(fn,'w') as f: f.create_dataset('pos',data=P.T); f.create_dataset('Cell',data=np.eye(3)*9); f.create_dataset('speed',data=V.T)
    r=MatReader(fn)
    assert isinstance(r.props['speed'],_H5Lazy) and isinstance(r.pos_array,_H5Lazy)
    eager=read_mat(fn); t=av.load(fn)
    assert t.n_frames==F
    for k in (0,3,5):
        assert np.array_equal(t.get(k).positions,P[:,:,k]) and np.array_equal(t.get(k)['speed'],V[:,k])
        assert np.array_equal(eager[k].positions,P[:,:,k])
