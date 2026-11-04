function a=modargs(mods)
% {'-m', spec, ...} from one modifier string or a cell array of them.
if ischar(mods)||isstring(mods), mods=cellstr(mods); end
a=reshape([repmat({'-m'},1,numel(mods));mods(:)'],1,[]);
end
