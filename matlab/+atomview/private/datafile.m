function [f,c]=datafile(data,keep)
% File name for data: a struct or matrix is saved to a temporary .mat (removed with c unless keep).
if nargin<2, keep=false; end; c=[];
if ischar(data)||isstring(data), f=char(data); return; end
if isnumeric(data), data=struct('Position',data); end
f=[tempname '.mat'];
w=whos('data'); v='-v7'; if w.bytes>1.8e9, v='-v7.3'; end
save(f,'-struct','data',v);
if ~keep, c=onCleanup(@() delete(f)); end
end
