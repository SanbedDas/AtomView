import csv
import os
import subprocess
import sys

import numpy as np
import pytest
from scipy.io import savemat

import atomview as av
from atomview import batch
from atomview.io import read_mat

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def fcc(n=4, a=3.6):
    basis = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    g = np.array([[i, j, k] for i in range(n) for j in range(n) for k in range(n)])
    return ((g[:, None, :] + basis[None]).reshape(-1, 3)) * a, n * a


def write_xyz(path, pos, L, frames=1, element="Cu"):
    with open(path, "w") as fh:
        for f in range(frames):
            fh.write(f'{len(pos)}\nLattice="{L} 0 0 0 {L} 0 0 0 {L}" Properties=species:S:1:pos:R:3 step={f}\n')
            for p in pos + 0.01 * f:
                fh.write(f"{element} {p[0]:.4f} {p[1]:.4f} {p[2]:.4f}\n")


def write_mat73(path, pos, L):
    import h5py
    with h5py.File(path, "w") as f:  # MATLAB -v7.3 layout: column-major, so stored transposed
        f.create_dataset("pos", data=pos.T)
        f.create_dataset("box", data=np.diag([L] * 3).T)
        f.create_dataset("type", data=np.ones((1, len(pos))))


@pytest.fixture
def tree(tmp_path):
    pos, L = fcc()
    (tmp_path / "runA").mkdir()
    (tmp_path / "runA" / "sub").mkdir()
    (tmp_path / "runB").mkdir()
    (tmp_path / "empty").mkdir()
    write_xyz(tmp_path / "runA" / "step_9.xyz", pos, L)
    write_xyz(tmp_path / "runA" / "step_10.xyz", pos, L, frames=3)  # natural order: 10 is last
    savemat(tmp_path / "runA" / "sub" / "final.mat",
            {"Position": pos, "Cell": np.diag([L] * 3), "type": np.arange(len(pos)) % 2 + 1,
             "energy": np.linspace(-3, -2, len(pos))[:, None]})
    savemat(tmp_path / "runB" / "a.mat", {"data": {"xyz": pos.T, "BoxSize": [L, L, L],
                                                    "species": np.array(["Cu"] * len(pos), dtype=object)}})
    write_mat73(tmp_path / "runB" / "b.mat", pos, L)
    return tmp_path


def test_latest_file_natural_and_mtime(tree):
    assert os.path.basename(batch.latest_file(str(tree / "runA"))) == "step_10.xyz"
    os.utime(tree / "runA" / "step_9.xyz", (2e9, 2e9))
    assert os.path.basename(batch.latest_file(str(tree / "runA"), pick="mtime")) == "step_9.xyz"
    assert batch.latest_file(str(tree / "empty")) is None


def test_mat_variants(tree):
    f = read_mat(str(tree / "runA" / "sub" / "final.mat"))[0]
    assert f.n == 256 and set(np.unique(f["type"])) == {1, 2} and "energy" in f.properties
    g = read_mat(str(tree / "runB" / "a.mat"))[0]  # struct + transposed positions + names
    assert g.n == 256 and g.type_names == {1: "Cu"} and np.allclose(np.diag(g.cell), 14.4)
    h = av.load(str(tree / "runB" / "b.mat")).get(0)  # v7.3 / HDF5
    assert h.n == 256 and np.allclose(np.diag(h.cell), 14.4)


def test_mat_trajectory(tmp_path):
    pos, L = fcc(2)
    traj = np.stack([pos + 0.1 * k for k in range(4)], axis=2)  # N x 3 x F
    savemat(tmp_path / "t.mat", {"Position": traj, "Cell": np.eye(3) * L, "Timestep": [0, 10, 20, 30]})
    t = av.load(str(tmp_path / "t.mat"))
    assert t.n_frames == 4 and t.get(3).attributes["Timestep"] == 30


def test_batch_cli(tree, tmp_path):
    out = tmp_path / "images"
    r = subprocess.run([sys.executable, "-m", "atomview", "batch", str(tree), "-o", str(out), "--size", "320x240"],
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout + r.stderr
    names = sorted(os.listdir(out))
    assert names == ["batch_summary.csv", "runA__step_10.png", "runA__sub__final.png", "runB__b.png"], names
    rows = list(csv.DictReader(open(out / "batch_summary.csv")))
    step10 = next(r for r in rows if r["file"].endswith("step_10.xyz"))
    assert step10["frame"] == "2"  # last frame of the multi-frame file


def test_batch_next_to_files_with_session(tree):
    p = av.Pipeline(av.load(str(tree / "runA" / "step_9.xyz")))
    p.add(av.modifiers.CommonNeighborAnalysis())
    session = tree / "s.json"
    p.save(session)
    res = batch.run(str(tree), session=str(session), size=(200, 150), progress=lambda m: None)
    assert all(r["status"] == "ok" for r in res), res
    assert os.path.exists(tree / "runB" / "b.png") and os.path.exists(tree / "runA" / "sub" / "final.png")


def test_frame_spec():
    fs = batch.frame_spec
    assert fs("last", 10) == [9] and fs("first", 10) == [0] and fs("all", 4) == [0, 1, 2, 3]
    assert fs("::5", 12) == [0, 5, 10] == fs("every5", 12) == fs("every 5", 12)
    assert fs("-3:", 10) == [7, 8, 9] and fs("0,5,-1,99", 10) == [0, 5, 9] and fs("-1", 10) == [9] and fs("50", 10) == [9]
    assert batch.single("last") and batch.single("-1") and not batch.single("::5") and not batch.single("0,5")


def snapshots(root):
    pos, L = fcc(2)
    (root / "a").mkdir(); (root / "b").mkdir()
    for k in range(12):
        write_xyz(root / "a" / f"snap_{k}.xyz", pos + 0.1 * k, L)  # one snapshot per file
    write_xyz(root / "b" / "traj.xyz", pos, L, frames=7)  # one file, seven frames


def test_batch_all_every_and_parallel(tmp_path):
    snapshots(tmp_path)
    r = batch.run(str(tmp_path), out_dir=str(tmp_path / "o1"), files="all", frames="every5", size=(160, 120), progress=lambda m: None)
    names = sorted(os.path.basename(x["image"]) for x in r)
    assert names == ["a__snap_0.png", "a__snap_10.png", "a__snap_5.png", "b__traj_00000.png", "b__traj_00005.png"]
    assert all(x["status"] == "ok" for x in r)
    r = batch.run(str(tmp_path), out_dir=str(tmp_path / "o2"), frames="all", size=(160, 120), progress=lambda m: None)
    assert len(r) == 8  # last file of a (1 frame) + all 7 frames of b
    r2 = batch.run(str(tmp_path), out_dir=str(tmp_path / "o3"), frames="all", size=(160, 120), jobs=2, progress=lambda m: None)
    assert [x["image"].replace("o3", "o2") for x in r2] == [x["image"] for x in r]  # same images, same order
    import imageio.v2 as iio
    d = [np.abs(iio.imread(x["image"]).astype(float) - iio.imread(y["image"])) for x, y in zip(r, r2)]
    assert max(e.mean() for e in d) < 0.5, [(e.mean(), e.max()) for e in d]  # GPU drivers may differ by a few pixels between processes
    r3 = batch.run(str(tmp_path), out_dir=str(tmp_path / "o2"), frames="all", size=(160, 120), skip_existing=True, progress=lambda m: None)
    assert all(x["status"].startswith("skipped") for x in r3) and len(r3) == 8
