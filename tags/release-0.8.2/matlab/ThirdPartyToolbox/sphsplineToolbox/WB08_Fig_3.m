% $Id: WB08_Fig_3.m,v 1.1.1.1 2008/05/09 21:34:52 myself Exp $
%
% This prepares data for Figure 3 which is created with the corresponding
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


%Then WB's solution (p = 38.9)
% p = 38.9;
p = 10;
Z = sphsplinet (loni, lati, zi, X, Y, p);
subplot (2,1,2)
contour (X, Y, Z)
A = [X(:) Y(:) Z(:)];
save Fig_2_p5.d A -ascii -tabs

% Validation

load mag_validate_1990.d
X = mag_validate_1990(:,1);
Y = mag_validate_1990(:,2);
Z = mag_validate_1990(:,3);
% First Parker's solution (p = 0)
Z0 = sphsplinet (loni, lati, zi, X, Y);
disp ('Parker differences:')
d0 = Z - Z0
rms(d0)
%Then WB's solution
Rmin = 1e100;
p_list = 36:0.1:42;
for p = p_list
    Zp = sphsplinet (loni, lati, zi, X, Y, p);
    dp = Z - Zp;
    R = rms (dp);
    if (R < Rmin)
        Rmin = R;
        pmin = p;
    end
end
disp ('WB08 differences:')
Rmin
pmin

