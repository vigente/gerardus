% $Id: WB08_Fig_4.m,v 1.1.1.1 2008/05/09 21:34:52 myself Exp $
%
% This prepares data for Figure 4 which is created with the corresponding
% shell script using GMT.  It also makes a plot in Matlab/Octave.
%
% Wessel, P. and J. M. Becker, 2008, Interpolation using a
%  generalized Green's function for a spherical surface spline
%  in tension, Geophys. J. Int., doi:10.1111/j.1365-246X.2008.03829.x
%
% Do the Mars example.

load mars370.in
loni = mars370(:,1);
lati = mars370(:,2);
zi   = mars370(:,3);
d=1;
% Set global 1x1 grid output coordinates

[X Y] = meshgrid (0:d:360, -90:d:90);

%Parker solution (p = 0)
Z = sphsplinet (loni, lati, zi, X, Y);
subplot (2,1,2)
contour (X, Y, Z)
A = [X(:) Y(:) Z(:)];
save Fig_4a.d A -ascii -tabs


%Then WBs solution (p = 20)
p = 20;
Z = sphsplinet (loni, lati, zi, X, Y, p);
subplot (2,1,2)
contour (X, Y, Z)
A = [X(:) Y(:) Z(:)];
save Fig_4b.d A -ascii -tabs

