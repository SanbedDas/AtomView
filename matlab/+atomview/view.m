function view(data)
% Open data (struct, matrix or file) in the AtomView window without blocking MATLAB.
[f,c]=datafile(data,true); %#ok<ASGLU>
avrun({f},'',true);
end
