% Checks the AtomView MATLAB toolbox end to end: run from the repository root with
%   matlab -batch "addpath('matlab'); run('matlab/tests/atomview_test.m')"
root=fileparts(fileparts(fileparts(mfilename('fullpath')))); d=tempname; mkdir(d); ok=@(c,m) assert(c,m);
fprintf('cli: %s\n',atomview.setup());
rng(1); N=400; F=30; L=40;
P=cumsum(cat(3,rand(N,2)*L,0.3*randn(N,2,F-1)),3);  % 2D random walk, N x 2 x F
S=struct('Position',[P,zeros(N,1,F)],'Cell',diag([L L 1]),'type',1+(rand(N,1)>0.5),'id',(1:N)','Timestep',(0:F-1)*10);
S.Position(:,1:2,:)=mod(S.Position(:,1:2,:),L);  % stored wrapped
% load round trip of a MATLAB struct
R=atomview.load(S); ok(isequal(size(R.Position),[N 3 F]) && max(abs(R.Position(:)-S.Position(:)))==0,'round trip');
fprintf('PASS load(struct) round trip, %d x %d x %d\n',size(R.Position));
% dynamics vs a direct MSD over all time origins, unwrapped walk
T=atomview.dynamics(S,1.5);
msd=zeros(F,1); for m=1:F-1, dr=P(:,:,1+m:F)-P(:,:,1:F-m); msd(m+1)=mean(sum(dr.^2,2),'all'); end
ok(max(abs(T.('MSD(t)')-msd))<1e-9*max(msd),'msd'); ok(T.time(2)==10,'time axis from Timestep');
fprintf('PASS dynamics: MSD matches direct MATLAB computation (max rel err %.1e), columns: %s\n',max(abs(T.('MSD(t)')-msd))/max(msd),strjoin(T.Properties.VariableNames,', '));
% render a struct, and a file with modifiers and an overlay
atomview.render(S,fullfile(d,'s.png'),{'ColorCoding(property=''type'')'},[500 400],[],'--frame 29');
atomview.render(fullfile(root,'examples','data','cu_void.dump'),fullfile(d,'ptm.png'),{'PolyhedralTemplateMatching()','Slice(width=8)'},[800 600],'-y','--overlay "t = {Timestep}"');
ok(exist(fullfile(d,'s.png'),'file')&&exist(fullfile(d,'ptm.png'),'file'),'render'); fprintf('PASS render struct and file -> %s\n',d);
% analysis results back into MATLAB
C=atomview.load(fullfile(root,'examples','data','cu_void.dump'),{'PolyhedralTemplateMatching()'},[0 20 5]);
ok(isequal(size(C.structure_type),[3921 4]) && mean(C.structure_type(:,1)==1)>0.9,'ptm');
fprintf('PASS load with PTM: FCC fraction per frame %s, orientation %s\n',mat2str(round(mean(C.structure_type==1),3)),mat2str(size(C.orientation)));
% batch over folders of .mat files written by MATLAB
for k=1:3, sub=fullfile(d,'runs',sprintf('run%d',k)); mkdir(sub); for t=[5 10], x=S; x.Position=S.Position(:,:,1:t); save(fullfile(sub,sprintf('step_%d.mat',t)),'-struct','x'); end; end
B=atomview.batch(fullfile(d,'runs'),fullfile(d,'imgs'),{},[300 200]);
ok(height(B)==3 && all(contains(B.file,'step_10')) && all(startsWith(B.status,'ok')),'batch');
fprintf('PASS batch: %d images, newest file per folder\n',height(B));
% movie and a Python script through the bundled interpreter
atomview.animate(S,fullfile(d,'walk.gif'),{'Trails(length=8)'},[300 300],[],'--every 3 --fps 8'); ok(exist(fullfile(d,'walk.gif'),'file')>0,'animate');
py=fullfile(d,'s.py'); fid=fopen(py,'w'); fprintf(fid,'import sys, atomview, numpy\nprint("py ok", atomview.__version__, sys.argv[1:])\n'); fclose(fid);
txt=atomview.script(py,'a1'); ok(contains(txt,'py ok'),'script'); fprintf('PASS animate and script: %s\n',strtrim(txt));
disp('ALL MATLAB TOOLBOX TESTS PASSED');
