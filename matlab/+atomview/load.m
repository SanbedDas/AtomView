function S=load(data,mods,frames)
% Read any file AtomView understands (LAMMPS, XYZ, POSCAR, XTC/DCD/GRO, .mat, ...) into a struct, after
% optional modifiers, e.g. atomview.load('dump.lmp',{'PolyhedralTemplateMatching()'}).
% frames: [start stop step] in 0-based frame numbers (stop exclusive), default all.
% Position is N x 3 x F; per-particle results N x F (or N x k x F); attr_<name> 1 x F; table_<name> [x y].
if nargin<2, mods={}; end; if nargin<3, frames=[]; end
[f,c]=datafile(data); %#ok<ASGLU>
o=[tempname '.mat']; k=onCleanup(@() delete(o)); %#ok<NASGU>
a=[{'convert',f,o},modargs(mods)];
if ~isempty(frames), fr=[frames(:)' 1]; a=[a,{'--frames',sprintf('%d:%d:%d',fr(1),fr(2),fr(3))}]; end
avrun(a,'');
S=builtin('load',o);
end
