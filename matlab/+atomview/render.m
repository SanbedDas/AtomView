function out=render(data,out,mods,sz,view,extra)
% Render one image. data: struct (Position, type, Cell, ... as in .mat files), N x 3 matrix or file name.
% mods: modifier strings, e.g. {'PolyhedralTemplateMatching()','Slice(width=6)'}; sz [w h]; view 'iso','+x',...
% extra: further CLI options, e.g. '--frame 5 --ao --overlay "t = {Timestep}"'.
if nargin<3, mods={}; end; if nargin<4||isempty(sz), sz=[1600 1200]; end
if nargin<5||isempty(view), view='iso'; end; if nargin<6, extra=''; end
[f,c]=datafile(data); %#ok<ASGLU>
avrun([{'render',f,'-o',out,'--size',sprintf('%dx%d',sz),['--view=' view]},modargs(mods)],extra);
end
