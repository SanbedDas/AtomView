function out=animate(data,out,mods,sz,view,extra)
% Render a movie (.mp4, .gif, .avi) of every frame. data: struct with Position N x 3 x F, matrix or trajectory file.
% extra e.g. '--fps 15 --every 2 --overlay "t = {Timestep*0.002:.1f} ps"'.
if nargin<3, mods={}; end; if nargin<4||isempty(sz), sz=[1280 960]; end
if nargin<5||isempty(view), view='iso'; end; if nargin<6, extra=''; end
[f,c]=datafile(data); %#ok<ASGLU>
avrun([{'animate',f,'-o',out,'--size',sprintf('%dx%d',sz),['--view=' view]},modargs(mods)],extra);
end
