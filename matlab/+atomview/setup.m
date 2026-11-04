function exe=setup(exe)
% atomview.setup('/path/to/atomview-cli') remembers the AtomView command; atomview.setup() shows the one in use.
if nargin, setpref('atomview','cli',exe); end
exe=cli();
end
