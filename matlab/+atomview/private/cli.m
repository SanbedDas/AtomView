function exe=cli()
% AtomView command: ATOMVIEW_CLI, atomview.setup, a checkout next to this toolbox, a standalone build, PATH.
exe=getenv('ATOMVIEW_CLI'); if ~isempty(exe), return; end
if ispref('atomview','cli'), exe=getpref('atomview','cli'); return; end
root=fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));  % private -> +atomview -> matlab -> repository
c={fullfile(root,'atomview-cli'),fullfile(root,'atomview-cli.exe'),fullfile(root,'AtomView.app','Contents','MacOS','atomview-cli'), ...
   fullfile(root,'.venv','bin','atomview'),fullfile(root,'.venv','Scripts','atomview.exe'), ...
   fullfile(root,'dist','AtomView','atomview-cli'),fullfile(root,'dist','AtomView','atomview-cli.exe'), ...
   fullfile(root,'build','dist','AtomView','atomview-cli')};
for i=1:numel(c), if exist(c{i},'file'), exe=c{i}; return; end; end
for n={'atomview-cli','atomview'}
    if ispc, [s,p]=system(['where ' n{1} ' 2>nul']); else, [s,p]=system(['which ' n{1} ' 2>/dev/null']); end
    if s==0, exe=strtrim(strtok(p,newline)); return; end
end
error('atomview:notFound','AtomView not found. Run atomview.setup(''/path/to/atomview-cli'') once.');
end
