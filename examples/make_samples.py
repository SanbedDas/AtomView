"""Generate small sample datasets in examples/data/.

    python examples/make_samples.py
"""
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "data")


def fcc(n, a):
    basis = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    g = np.array([[i, j, k] for i in range(n) for j in range(n) for k in range(n)])
    return ((g[:, None, :] + basis[None]).reshape(-1, 3)) * a


def bcc(n, a):
    basis = np.array([[0, 0, 0], [0.5, 0.5, 0.5]])
    g = np.array([[i, j, k] for i in range(n) for j in range(n) for k in range(n)])
    return ((g[:, None, :] + basis[None]).reshape(-1, 3)) * a


def write_dump(path, frames, L):
    with open(path, "w") as fh:
        for step, (pos, types, pe) in frames:
            n = len(pos)
            fh.write(f"ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF ATOMS\n{n}\n")
            fh.write(f"ITEM: BOX BOUNDS pp pp pp\n0 {L}\n0 {L}\n0 {L}\n")
            fh.write("ITEM: ATOMS id type x y z c_pe\n")
            for i in range(n):
                fh.write(f"{i + 1} {types[i]} {pos[i, 0]:.5f} {pos[i, 1]:.5f} {pos[i, 2]:.5f} {pe[i]:.5f}\n")


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(0)
    a, n = 3.615, 10
    L = a * n
    base = fcc(n, a)
    center = np.full(3, L / 2)
    keep = np.linalg.norm(base - center, axis=1) > 6.0  # spherical void
    base = base[keep]
    types = np.where(base[:, 0] < L / 2, 1, 2)
    frames = []
    for f in range(20):
        amp = 0.03 + 0.012 * f
        pos = (base + rng.normal(0, amp, base.shape)) % L
        r = np.linalg.norm(base - center, axis=1)
        pe = -3.54 + 1.2 * np.exp(-(r - 6.0) / 1.5) + rng.normal(0, 0.02 + 0.01 * f, len(base))
        frames.append((f * 1000, (pos, types, pe)))
    write_dump(os.path.join(OUT, "cu_void.dump"), frames, L)

    # A BCC iron crystal as extended XYZ with element names
    p = bcc(8, 2.866) + rng.normal(0, 0.04, (2 * 8 ** 3, 3))
    L2 = 2.866 * 8
    with open(os.path.join(OUT, "fe_bcc.xyz"), "w") as fh:
        fh.write(f"{len(p)}\nLattice=\"{L2} 0 0 0 {L2} 0 0 0 {L2}\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"\n")
        for q in p % L2:
            fh.write(f"Fe {q[0]:.5f} {q[1]:.5f} {q[2]:.5f}\n")
    print("Wrote", os.listdir(OUT))


if __name__ == "__main__":
    main()
