% $Id: WB08_Fig_5.m,v 1.1.1.1 2008/05/09 21:34:52 myself Exp $
%
% This prepares data for Figure 5 which is created with the corresponding
% shell script using GMT.  It also makes a plot in Matlab/Octave.
%
% Wessel, P. and J. M. Becker, 2008, Interpolation using a
%  generalized Green's function for a spherical surface spline
%  in tension, Geophys. J. Int., doi:10.1111/j.1365-246X.2008.03829.x
%
% Do the Moon example.

load lun2.xyz
loni = lun2(:,1);
lati = lun2(:,2);
zi   = lun2(:,3);
% Set global 15x15 arc min grid output coordinates

d=0.25;
[X Y] = meshgrid (0:d:360, -90:d:90);

%Parker solution (p = 0)
Z = sphsplinet (loni, lati, zi, X, Y);
subplot (2,1,2)
contour (X, Y, Z)
A = [X(:) Y(:) Z(:)];
save Fig_5a.d A -ascii -tabs


%Then WBs solution (p = 25)
p = 25;
Z = sphsplinet (loni, lati, zi, X, Y, p);
subplot (2,1,2)
contour (X, Y, Z)
A = [X(:) Y(:) Z(:)];
save Fig_5b.d A -ascii -tabs

