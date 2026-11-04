# Scripting AtomView from Python and MATLAB

Everything here is free: rendering, batch jobs, analysis and the Python interpreter that runs your scripts.

## Python

```python
import atomview as av
from atomview import modifiers as m
from atomview.render import RenderSettings

pipe = av.Pipeline(av.load("dump.lammpstrj"))            # or .mat, .xyz, .xtc + topology, POSCAR, ...
pipe.add(m.PolyhedralTemplateMatching())
pipe.add(m.Slice(normal=(0, 1, 0), width=8))
pipe.add(lambda frame, data: data.attributes.update(FCC=int((data["structure_type"] == 1).sum())))

style = RenderSettings(ambient_occlusion=True, overlay="t = {Timestep*0.002:.1f} ps   FCC = {FCC}")
av.render_image(pipe, "frame0.png", frame=0, size=(2400, 1800), settings=style, view="-y")
av.render_animation(pipe, "movie.mp4", settings=style, view="-y")

data = pipe.compute(10)                # Frame: data.positions, data["structure_type"], data.attributes, data.tables
av.write(data, "frame10.mat")          # MATLAB-readable

from atomview.trajectory import dynamics
res = dynamics(pipe, k=3.0)            # MSD, alpha2, Fs(k,t), VACF tables
```

**Without a Python installation** (standalone build): `atomview-cli script my_script.py arg1 arg2` runs the
script with atomview, NumPy, SciPy, freud and ASE available.

**From the GUI**: add a *Python script* modifier (`def modify(frame, data): ...`) to the pipeline.

## MATLAB

Add the toolbox once: `addpath('/path/to/AtomView/matlab')`. It finds AtomView next to itself (standalone
build or checkout); otherwise run `atomview.setup('/path/to/atomview-cli')` once. MATLAB's Python engine is not
needed. Frames are numbered from 0, as in the AtomView window.

```matlab
% your simulation data straight from the workspace
S.Position = pos;            % N x 3, N x 2 (2D) or N x 3 x F
S.Cell     = diag([Lx Ly Lz]);
S.type     = types;          % plus any per-particle arrays: S.energy, S.theta, S.velocity, ...

atomview.render(S, 'fig.png')                                         % image
atomview.render(S, 'fig.png', {'ColorCoding(property=''energy'')'}, [2400 1800], '+z', '--ao')
atomview.animate(S, 'movie.mp4', {'Trails(length=20)', 'Vectors(property=''theta'')'}, [], [], '--fps 15')
atomview.view(S)                                                      % open it in the AtomView window

% analysis results back into MATLAB
R = atomview.load('dump.lmp', {'PolyhedralTemplateMatching()'});     % R.structure_type is N x F, R.orientation N x 4 x F
T = atomview.dynamics(S, 7.0);                                        % table: time, MSD, alpha2, Fs(k,t), VACF
B = atomview.batch('/data/sims', '/data/figures');                    % newest file per folder -> image; summary table
B = atomview.batch('/data/sims', '/data/figures', {}, [], [], '--frames=every5 -j 0');  % every 5th frame, all cores
atomview.script('analysis.py', 'arg1');                               % run Python from MATLAB
```

Every function takes `(data, out, mods, sz, view, extra)` in that order where it applies: `data` is a struct,
an N x 3 matrix or any file AtomView reads; `mods` is one modifier string or a cell array of them; `[]` keeps a
default; `extra` passes any command-line option (`--frame 5`, `--overlay "..."`, `--transparent`, `--topology top.gro`).

The other direction also works: the *MATLAB / Octave script* modifier runs a `.m` file on every frame inside an
AtomView pipeline (one MATLAB session stays open, ~0.1 s per frame).

## Modifier strings

Any modifier class with keyword arguments, e.g. `ColorCoding(property='c_pe', colormap='plasma')`,
`ExpressionSelect(expression='type == 2 && z > 10')`, `CommonNeighborAnalysis()`, `AtomicStrain(cutoff=3.0)`,
`StructureFactor(k_max=12)`. Values must be literals; `atomview plugins` lists every modifier and its class name.
The same strings work on the command line: `atomview render dump.lmp -o fig.png -m "PolyhedralTemplateMatching()"`.
