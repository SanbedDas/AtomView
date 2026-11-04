"""Analyses over a whole trajectory (as opposed to per-frame modifiers)."""
from __future__ import annotations

import numpy as np

from .data import ID


def unwrapped_positions(pipe, frames, only_selected=False, progress=None):
    # (F, N, 3) positions of the particles present in every frame (matched by id), unwrapped across
    # periodic boundaries by summing minimum-image steps between frames (valid while steps < half the cell).
    # Returns (positions, first frame, number of particles dropped because they were missing in some frame).
    from .data import Frame
    snaps, f0 = [], None
    for k, fr in enumerate(frames):
        if progress and progress(k, len(frames)) is False:
            from .pipeline import Cancelled
            raise Cancelled()
        f = pipe.compute(fr)
        if f0 is None:
            f0 = f
            keep = f.selection if only_selected else np.ones(f.n, bool)
        ids = f.properties[ID] if ID in f.properties else None
        snaps.append((ids, f.positions, Frame(np.zeros((0, 3)), f.cell, pbc=f.pbc)))
    if snaps[0][0] is not None:
        common = set(snaps[0][0][keep].tolist())
        for ids, _, _ in snaps[1:]:
            common &= set(ids.tolist())
        want = np.array(sorted(common), dtype=np.int64)
        dropped = int(keep.sum()) - len(want)
        rows = [np.searchsorted(ids, want, sorter=np.argsort(ids)) for ids, _, _ in snaps]
        rows = [np.argsort(ids)[r] for (ids, _, _), r in zip(snaps, rows)]
    else:
        if any(len(p) != len(snaps[0][1]) for _, p, _ in snaps): raise ValueError("Particle count changes and there is no id property to match particles")
        rows, dropped = [np.nonzero(keep)[0]] * len(snaps), 0
    if not len(rows[0]): raise ValueError("No particle is present in every frame")
    out = np.empty((len(frames), len(rows[0]), 3))
    for k, ((_, p, cellf), r) in enumerate(zip(snaps, rows)):
        out[k] = p[r] if k == 0 else out[k - 1] + cellf.minimum_image(p[r] - snaps[k - 1][1][rows[k - 1]])
    return out, f0, dropped


def windowed_msd(r):
    # MSD(m) averaged over all time origins and particles, float64 FFT method; r is (F, N, d).
    F, N = r.shape[:2]
    D = np.vstack([(r ** 2).sum(-1), np.zeros((1, N))])
    Q, S1 = 2 * D.sum(0), np.empty((F, N))
    for m in range(F):
        Q = Q - D[m - 1] - D[F - m]
        S1[m] = Q / (F - m)
    X = np.fft.rfft(r, n=2 * F, axis=0)
    S2 = np.fft.irfft((X * X.conj()).real.sum(-1), n=2 * F, axis=0)[:F] / (F - np.arange(F))[:, None]
    return (S1 - 2 * S2).mean(1)


def _check_frames(frames):
    if len(frames) < 2: raise ValueError("Needs at least two frames")
    if np.ptp(np.diff(frames)) != 0: raise ValueError("Needs equally spaced frames")


def time_step(pipe, frames, dt=None):
    # Time between analysed frames: dt per source frame times the stride, else from the Timestep attribute, else the stride.
    stride = frames[1] - frames[0]
    if dt is not None:
        return dt * stride
    t0, t1 = pipe.compute(frames[0]).attributes.get("Timestep"), pipe.compute(frames[1]).attributes.get("Timestep")
    return float(t1 - t0) if isinstance(t0, (int, float)) and isinstance(t1, (int, float)) and t1 != t0 else float(stride)


def msd_table(pos, f0, dropped, step, fit_from=0.5):
    m = windowed_msd(pos)
    m[0] = 0.0
    dim = 2 if f0.attributes.get("Dimensions") == 2 else 3
    lag = np.arange(len(m)) * step
    i0 = max(1, int(len(m) * fit_from))
    slope = np.polyfit(lag[i0:], m[i0:], 1)[0] if len(m) - i0 >= 2 else np.nan
    return {"x": lag, "y": m, "xlabel": "time", "ylabel": "MSD", "kind": "line", "D": slope / (2 * dim), "dim": dim,
            "particles": pos.shape[1], "dropped": dropped, "frames": len(pos)}


def msd(pipe, frames=None, dt=None, only_selected=False, fit_from=0.5, progress=None):
    """Time-origin averaged MSD(t). dt = time per source frame; default from Timestep attributes, else 1.
    Returns a table dict with x (time), y (MSD), D (diffusion coefficient from a linear fit of
    the lags after `fit_from` of the range) and dim."""
    frames = list(range(pipe.n_frames)) if frames is None else list(frames)
    _check_frames(frames)
    pos, f0, dropped = unwrapped_positions(pipe, frames, only_selected, progress)
    return msd_table(pos, f0, dropped, time_step(pipe, frames, dt), fit_from)


def lag_stats(pos, k=None, max_origins=200):
    # Per lag m: <dr^2>, <dr^4> and (if k) the isotropic self-intermediate scattering function, averaged over
    # particles and up to max_origins evenly spaced time origins.
    from scipy.special import j0
    F, N, dim = pos.shape
    r2, r4, fs = np.zeros(F), np.zeros(F), np.ones(F)
    for m in range(1, F):
        o = np.arange(0, F - m, max(1, (F - m) // max_origins))
        dr = np.linalg.norm(pos[o + m] - pos[o], axis=-1)
        r2[m], r4[m] = (dr ** 2).mean(), (dr ** 4).mean()
        if k:
            fs[m] = (j0(k * dr) if dim == 2 else np.sinc(k * dr / np.pi)).mean()
    return r2, r4, fs


def autocorrelation(v):
    # <v(t0) . v(t0 + m)> over all time origins and particles (FFT); v is (F, N, d).
    F = len(v)
    X = np.fft.rfft(v, n=2 * F, axis=0)
    return (np.fft.irfft((X * X.conj()).real.sum(-1), n=2 * F, axis=0)[:F] / (F - np.arange(F))[:, None]).mean(1)


def dynamics(pipe, frames=None, dt=None, k=None, only_selected=False, progress=None, max_origins=200):
    """MSD(t), non-Gaussian parameter alpha2(t), self-intermediate scattering Fs(k,t) (if k) and the normalized
    velocity autocorrelation, as tables keyed by name. Velocities come from a 'velocity' property in every frame,
    otherwise from finite differences of the unwrapped positions."""
    frames = list(range(pipe.n_frames)) if frames is None else list(frames)
    _check_frames(frames)
    pos, f0, dropped = unwrapped_positions(pipe, frames, only_selected, progress)
    step = time_step(pipe, frames, dt)
    r = msd_table(pos, f0, dropped, step)
    dim = r["dim"]
    pos = pos[..., :dim]
    t = r["x"]
    r2, r4, fs = lag_stats(pos, k, max_origins)
    with np.errstate(invalid="ignore", divide="ignore"):
        a2 = dim * r4 / ((dim + 2) * r2 ** 2) - 1
    a2[0] = np.nan
    out = {"MSD(t)": r,
           "alpha2(t)": {"x": t, "y": a2, "xlabel": "time", "ylabel": "alpha2", "kind": "line"}}
    if k:
        out[f"Fs(k={k:g},t)"] = {"x": t, "y": fs, "xlabel": "time", "ylabel": "Fs(k,t)", "kind": "line"}
    vel = [pipe.compute(fr).get("velocity") for fr in frames]
    if all(v is not None and len(v) == pos.shape[1] for v in vel) and f0.n == pos.shape[1]:
        v, src = np.stack(vel)[..., :dim].astype(float), "velocity property"
    else:
        v, src = np.diff(pos, axis=0) / step, "finite differences"
    c = autocorrelation(v)
    out["VACF(t)"] = {"x": np.arange(len(c)) * step, "y": c / c[0] if c[0] else c, "xlabel": "time",
                      "ylabel": "C_v(t)/C_v(0)", "kind": "line", "source": src}
    return out
