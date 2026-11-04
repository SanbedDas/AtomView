% MATLAB / Octave script modifier example: load this file in a
% "MATLAB / Octave script" modifier. Runs once per frame.
center = Origin + 0.5 * sum(Cell, 1);
d = Position - center;
distance_from_center = sqrt(sum(d.^2, 2));     % new N x 1 property
selection = distance_from_center < 9;          % select atoms near the void
NearVoid = sum(selection);                     % scalar -> global attribute
