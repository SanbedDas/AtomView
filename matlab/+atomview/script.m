function txt=script(py,varargin)
% Run a Python script with atomview, numpy, scipy, freud and ase available (also from a standalone build).
txt=avrun([{'script',py},varargin],'');
end
