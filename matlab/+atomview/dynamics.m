function T=dynamics(data,k,every,dt,mods)
% MSD, non-Gaussian parameter alpha2, Fs(k,t) (if k given) and velocity autocorrelation as a table.
if nargin<2, k=[]; end; if nargin<3||isempty(every), every=1; end; if nargin<4, dt=[]; end; if nargin<5, mods={}; end
[f,c]=datafile(data); %#ok<ASGLU>
o=[tempname '.csv']; x=onCleanup(@() delete(o)); %#ok<NASGU>
a=[{'dynamics',f,'-o',o,'--every',num2str(every)},modargs(mods)];
if ~isempty(k), a=[a,{'--k',num2str(k,17)}]; end
if ~isempty(dt), a=[a,{'--dt',num2str(dt,17)}]; end
avrun(a,'');
T=readtable(o,'VariableNamingRule','preserve');
end
