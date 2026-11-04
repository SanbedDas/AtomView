function txt=avrun(args,extra,background)
% Calls the AtomView command with quoted arguments; errors with its output if it fails.
if nargin<3, background=false; end
q=cellfun(@quote,[{cli()},args],'UniformOutput',false); cmd=strjoin(q,' ');
if ~isempty(extra), cmd=[cmd ' ' extra]; end
if ispc
    if background, cmd=['start "" ' cmd]; end
else
    cmd=['env -u LD_LIBRARY_PATH -u QT_PLUGIN_PATH -u PYTHONHOME -u PYTHONPATH ' cmd];  % MATLAB sets these for its own libraries
    if background, cmd=[cmd ' > /dev/null 2>&1 &']; end
end
[s,txt]=system(cmd);
txt=strjoin(regexp(txt,'^(?!.*[Ff]ontconfig)(?!Matplotlib is building)(?!\s+/).*$','match','lineanchors'),newline);
if s~=0, error('atomview:failed','%s\n%s',cmd,txt); end
end

function a=quote(a)
a=char(string(a));
if ispc, a=['"' strrep(a,'"','\"') '"']; else, a=['''' strrep(a,'''','''\''''') '''']; end
end
