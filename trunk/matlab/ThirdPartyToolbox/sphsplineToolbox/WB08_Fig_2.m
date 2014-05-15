% $Id: WB08_Fig_2.m,v 1.1.1.1 2008/05/09 21:34:52 myself Exp $
%
% This prepares data for Figure 2 which is created with the corresponding
% shell script using GMT.  It also makes a plot in Matlab/Octave.
%
% Wessel, P. and J. M. Becker, 2008, Interpolation using a
%  generalized Green's function for a spherical surface spline
%  in tension, Geophys. J. Int., doi:10.1111/j.1365-246X.2008.03829.x
%
% Replicate Parker and then find the tension that minimizes the misfit at
% his 8 extra validation stations.

load mag_obs_1990.d
loni = mag_obs_1990(:,1);
lati = mag_obs_1990(:,2);
zi   = mag_obs_1990(:,3);
d=1;
% Set global 1x1 grid output coordinates

[X Y] = meshgrid (0:d:360, 0:d:90);

% First Parker's solution (p = 0)
Z = sphsplinet (loni, lati, zi, X, Y);

figure(1); clf
subplot (2,1,1)
contour (X, Y, Z)
drawnow
A = [X(:) Y(:) Z(:)];
save Fig_2_p0.d A -ascii -tabs

%Then used the wrong Oslo longitude to recreate Parker's figure
k = find (loni == 10.45)
loni(k) = 104.5;
Z = sphsplinet (loni, lati, zi, X, Y);
subplot (2,1,2)
contour (X, Y, Z)
A = [X(:) Y(:) Z(:)];
save Fig_2_orig.d A -ascii -tabs
