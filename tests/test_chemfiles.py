import numpy as np
import chemfiles as cf
import atomview as av

def write(d,N=30,F=4):
    rng=np.random.default_rng(0); P=rng.uniform(1,20,(F,N,3)); c=cf.UnitCell([22,24,26],[90,85,75])
    with cf.Trajectory(str(d/'w.gro'),'w') as tg, cf.Trajectory(str(d/'w.xtc'),'w') as tx, cf.Trajectory(str(d/'w.dcd'),'w') as td:
        for k in range(F):
            fr=cf.Frame(); fr.cell=c; fr.step=k*10
            for i in range(N): fr.add_atom(cf.Atom('OW' if i%3==0 else 'HW'),P[k,i])
            if k==0: tg.write(fr)
            tx.write(fr); td.write(fr)
    return P,np.array(c.matrix).T

def test_xtc_dcd_gro_with_sibling_topology(tmp_path):
    P,cell=write(tmp_path)
    for fn,tol in (('w.xtc',0.006),('w.dcd',1e-5),('w.gro',0.006)):
        t=av.load(str(tmp_path/fn)); k=t.n_frames-1; f=t.get(k)
        assert t.format=='chemfiles' and np.abs(f.positions-P[k]).max()<tol and np.allclose(f.cell,cell,atol=0.01)
        assert f.type_names=={1:'OW',2:'HW'} and f.pbc==(True,True,True)
    assert av.load(str(tmp_path/'w.xtc')).n_frames==4 and av.load(str(tmp_path/'w.xtc')).get(3).attributes['Timestep']==30

def test_explicit_topology_and_pipeline(tmp_path):
    P,_=write(tmp_path); (tmp_path/'w.gro').rename(tmp_path/'top.gro')
    assert av.load(str(tmp_path/'w.xtc')).get(0).type_names!={1:'OW',2:'HW'}  # no sibling topology now
    t=av.load(str(tmp_path/'w.xtc'),topology=str(tmp_path/'top.gro')); assert t.get(0).type_names=={1:'OW',2:'HW'}
    p=av.Pipeline(t); p.add(av.modifiers.SelectType(types='OW')); assert p.compute(2).selection.sum()==10

def test_pdb_bonds_and_residues(tmp_path):
    fr=cf.Frame()
    for i,(n,x) in enumerate([('C',(0,0,0)),('O',(1.2,0,0)),('N',(2.4,0,0))]): fr.add_atom(cf.Atom(n),x)
    r=cf.Residue('ALA',1)
    for i in range(2): r.atoms.append(i)
    fr.add_residue(r); fr.add_bond(0,1)
    with cf.Trajectory(str(tmp_path/'m.pdb'),'w') as t: t.write(fr)
    f=av.load(str(tmp_path/'m.pdb')).get(0)
    assert f.n==3 and len(f.bonds)>=1 and f.pbc==(False,False,False) and 'residue' in f.properties
