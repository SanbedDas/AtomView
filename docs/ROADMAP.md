# AtomView: assessment, gaps, architecture and roadmap

Status as of 2026-11-03, commit `daa7e86` plus the work listed under "Done". Every number below was measured on an
Apple-silicon Mac (Python 3.14, VTK 9.7.1) with `benchmarks/bench.py`. Nothing here is a projection unless it is
labeled as a target.

## 1. Repository assessment

| Area | State |
|---|---|
| Size | ~4.7k lines of Python in 16 modules, 21 tests (core, batch, 2D), GUI smoke test, CI on Ubuntu/Windows/macOS (green) |
| Data model | `Frame` (positions, cell, pbc, named per-particle arrays, attributes, bonds, tables) - simple and sufficient |
| IO | LAMMPS dump/data, (ext)XYZ, POSCAR, `.mat` v5/v7/v7.3, ASE fallback. Lazy byte-offset index for dump/XYZ |
| Pipeline | Source + ordered modifiers, per-frame cache, JSON sessions, errors contained per modifier |
| Analysis | a-CNA (verified on FCC/HCP/BCC), CSP, RDF, clusters, displacement/MSD (per frame), Steinhardt, Voronoi, bonds, hexatic (2D) |
| Scripting | Python modifier, MATLAB/Octave modifier (subprocess per frame), Python API, CLI, batch |
| Rendering | PyVista sphere *meshes*, legends/colorbars, AO, SVG/PDF, MP4/GIF |

Measured baseline (offscreen 1600x1200 render of an FCC crystal):

| N | load frame | a-CNA | render | peak memory |
|---|---|---|---|---|
| 11 k | 0.003 s | 0.06 s | 1.21 s | 0.75 GB |
| 98 k | 0.03 s | 0.53 s | 1.38 s | 3.7 GB |
| 1.0 M | - | - | **no image: process died** | **15.5 GB** |

Problems found during assessment:

1. **Rendering does not scale.** Each atom becomes ~300 triangles of real geometry; 1M atoms exhausts memory.
2. **`.mat` trajectories load fully into memory** (every frame parsed up front), so long v7.3 trajectories cost
   RAM proportional to the whole file.
3. **No 2D support.** N x 2 positions were rejected; 2D RDF/neighbor normalization was wrong. *(Fixed in increment 1.)*
4. **Trajectory-level analysis is missing.** MSD exists only as "displacement from one reference frame"; there is
   no MSD(t) curve, the basic observable for colloids and active matter.
5. **MATLAB modifier costs ~6 s/frame** (fresh process). One cold start stalled >300 s with no feedback.
6. **No vector glyphs** (velocities, orientations of active particles) and **no plugin mechanism** for
   third-party modifiers.

## 2. Feature gap analysis

| Capability | AtomView | Priority |
|---|---|---|
| 2D systems (colloids, active matter) | done (increment 1) | P0 |
| Large systems (1M+ particles) interactive and exported | missing | P0 |
| Lazy `.mat` v7.3 trajectories | missing | P0 |
| MSD(t), time series over a trajectory | missing | P0 |
| Vector/arrow glyphs (velocity, orientation) | missing | P1 |
| Plugin modifiers (entry points) | missing | P1 |
| Persistent MATLAB/Octave engine | **done**: file-based job loop (no matlab.engine, which needs Python <= 3.13); 0.09 s/frame vs 8.3 s | P1 |
| a-CNA memory at 1M particles | **done**: CNA+CSP 8.97 GB -> 0.66 GB, 13.6 s -> 8.5 s, results identical | P1 |
| Progress and cancel for long computations | **done** (pipeline, CNA, MATLAB process, MSD); a running Python function finishes its frame first | P1 |
| Provenance (pipeline JSON embedded in exported images) | missing | P1 |
| Polyhedral template matching | **done**: bundled MIT reference library as a C++ extension; 1M atoms in 4.4 s | P2 |
| Dislocation analysis (DXA), surface mesh | missing | P3 (large; re-implement from papers only) |
| Molecular representations (cartoon/ribbons, PDB topology) | missing | P3 |
| Trajectory smoothing, unwrapping across frames | missing | P2 |
| Voronoi, clusters, CNA, CSP, RDF, Steinhardt, scripting, batch, headless | present | - |

## 3. Architecture (current to target)

```
sources (io.Trajectory, read_mat, ASE)  ->  Frame  ->  Pipeline[modifiers]  ->  render / GUI / CLI / batch
                                                            |
                                     trajectory analyses (multi-frame, new)  ->  tables / plots
```

Changes, each kept backward compatible:

- **Render backends.** `RenderSettings.render_mode`: `mesh` (current, best quality, small N), `imposter`
  (GPU point-sprite spheres via `vtkPointGaussianMapper`, constant memory per particle), `auto` (switch at a
  particle threshold). Existing sessions keep working.
- **Lazy sources.** `.mat` v7.3 frames read by HDF5 hyperslab on demand; v5 files (which SciPy cannot
  slice) stay eager.
- **Trajectory analyses.** Functions over `(pipeline, frames)` that return tables (MSD(t) first), callable
  from the API, CLI and GUI Plots tab.
- **Plugins.** Modifiers registered through the `atomview.modifiers` entry-point group, plus a user plugin
  folder, all through the existing `REGISTRY`.

## 4. Dependency and license risks

| Dependency | License | Risk / mitigation |
|---|---|---|
| VTK, PyVista, pyvistaqt | BSD-3 / MIT | Wheels needed for each new Python; pin in CI. Headless Linux depends on OSMesa (`libosmesa6`) |
| PySide6 (Qt) | LGPL-3.0 | Fine for dynamic linking from Python. A frozen redistributable bundle must let users swap Qt libraries |
| ASE | LGPL-2.1+ | Imported library, unmodified: compatible with MIT application code |
| freud | BSD-3 | Primary neighbor/order engine; API changed between majors, so pin `<4` |
| h5py, NumPy, SciPy, Matplotlib, imageio | BSD / PSF-style | Low |
| imageio-ffmpeg | BSD-2 wrapper, bundled ffmpeg binary | Check the bundled ffmpeg build's license before redistributing binaries |
| MATLAB | proprietary, optional | Never required; Octave (GPL, separate process) as the free path |

## 5. Prioritized roadmap with acceptance criteria

| # | Increment | Acceptance criteria (measurable) | Status |
|---|---|---|---|
| 1 | 2D data | N x 2 / N x 2 x F `.mat` and flat-z files load as 2D; periodic triangular lattice gives coordination 6 and abs(psi6)=1 to 1e-4; ideal 2D gas g(r) within 3% of 1 | **done** |
| 2 | Lazy `.mat` v7.3 | Opening an F-frame v7.3 trajectory reads only the requested frame; frame-k positions equal the eager reader exactly; memory to open a 2 GB file < 200 MB | **done**: 458 MB file, 932 MB -> 27 MB extra RSS, identical frames on 7 layouts |
| 3 | Imposter rendering | 1M-particle image renders with peak memory < 4 GB; `auto` keeps meshes below the threshold so small-system images are unchanged | **done**: 1M atoms render in 1.15 s using +0.16 GB (was: crash at 15.5 GB); 98k: 3.56 GB -> 0.81 GB; image coverage matches mesh within 10% |
| 4 | MSD(t) | Ballistic test (constant velocity) gives MSD = v^2 t^2 exactly; 2D random walk gives slope 4D within 5% | **done**: ballistic through periodic wrapping to 1e-9; FFT = direct sum to 1e-10; 2D D within 5%; freud's float32 MSD rejected (3e-5 error) |
| 5 | Vectors + plugins | Arrow glyphs from any N x 3 / N x 2 property; a modifier in an installed package's entry point appears in the GUI menu with no AtomView code change | **done**: arrows from N x 3, N x 2 or angle; entry-point and folder plugins load, a broken plugin is reported, not fatal |
| 6 | Persistent MATLAB engine | Per-frame cost after the first < 1 s; same outputs as the subprocess path | **done**: 0.09 s/frame, identical outputs, isolated workspaces, cancel kills and restarts |
| 7 | Provenance | Exported PNG carries the pipeline JSON and versions in its metadata; `atomview info fig.png` prints it | |
| 8 | PTM | Matches a-CNA on perfect FCC/HCP/BCC; RMSD output per particle | |

## Done after the first roadmap (2026-11-03)

| Item | Validation |
|---|---|
| Per-frame text overlays | `{expr:fmt}` over attributes; frame 12 of the sample renders "t = 24.0 ps" |
| chemfiles formats (XTC, TRR, DCD, NetCDF, GRO, PDB, ...) | GRO/XTC round trip within 0.005 A, DCD exact, triclinic cell, sibling topology |
| Progress and cancel | pipeline, CNA chunks, MATLAB process killed within seconds, MSD; GUI busy bar tested |
| CNA/CSP memory | 1M atoms: 8.97 GB -> 0.66 GB, identical results on 24,550 values |
| Atomic strain | simple shear and expansion exact to 1e-12; D2min local to a moved atom; 2D |
| S(k), polar/nematic, g6(r), trails | gas S=1, FCC (111) peak; exact order limits; g6=1 on a lattice |
| alpha2, Fs(k,t), VACF | same-speed ballistic: alpha2 = d/(d+2)-1, Fs = sinc / J0 exact to 1e-9 |
| PTM | all 7 lattices exact (RMSD ~1e-8, lattice constants exact); orientation recovered under an arbitrary rotation; Cu3Au L1_2 ordering; superset of CNA-FCC on the thermal sample |
