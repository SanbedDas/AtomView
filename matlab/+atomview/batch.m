function r=batch(folder,outdir,mods,sz,view,extra)
% Images of the snapshots in every sub-folder of folder (default: last frame of the last file); returns the summary table.
% extra e.g. '--frames=every5 --files all -j 0 --skip-existing --pattern "*.mat,*.dump"'.
if nargin<2, outdir=''; end; if nargin<3, mods={}; end; if nargin<4||isempty(sz), sz=[2000 1500]; end
if nargin<5||isempty(view), view='iso'; end; if nargin<6, extra=''; end
a={'batch',folder,'--size',sprintf('%dx%d',sz),['--view=' view]};
if ~isempty(outdir), a=[a,{'-o',outdir}]; end
avrun([a,modargs(mods)],extra);
s=fullfile(folder,'batch_summary.csv'); if ~isempty(outdir), s=fullfile(outdir,'batch_summary.csv'); end
r=readtable(s,'Delimiter',',','TextType','string');
end
