% $Id: WB08_Fig_1.m,v 1.1.1.1 2008/05/09 21:34:52 myself Exp $
%
% This prepares data for Figure 1 which is created with the corresponding
% shell script using GMT.  It also makes a plot in Matlab/Octave.
%
% Wessel, P. and J. M. Becker, 2008, Interpolation using a
%  generalized Green's function for a spherical surface spline
%  in tension, Geophys. J. Int., doi:10.1111/j.1365-246X.2008.03829.x
%
% This script plots the Green function in top panel and the gradient in the
% bottom panel, for tensions p = 0, 1, 5, 25, and 100.  We do this for all
% angles between 0 and 180 for 1001 equidistant points.
% The results are also written to files Fig_1_[g]p###.d (g for gradient)
% for making Figure 1 in WB08_Fig_1.sh

disp ('Do surface')
theta = linspace (0, 180, 1001)';
x = cosd (theta);
y = SSSTc (x, 0);
figure(2)
clf
subplot (2,1,1)
plot (theta, y)
A = [theta y];
save Fig_1_p0.d A -ascii -tabs
hold on
y = SSSTc (x, 1);
plot (theta, y)
A = [theta y];
save Fig_1_p1.d A -ascii -tabs
y = SSSTc (x, 5);
plot (theta, y)
A = [theta y];
save Fig_1_p5.d A -ascii -tabs
y = SSSTc (x, 25);
plot (theta, y)
A = [theta y];
save Fig_1_p25.d A -ascii -tabs
y = SSSTc (x, 100);
plot (theta, y)
A = [theta y];
save Fig_1_p100.d A -ascii -tabs

disp ('do gradient')
subplot (2,1,2)
y = real(SSSTgradc (x, 0));
k = find (y ~= 0);
y(k) = -y(k)/min(y);
plot (theta, y)
A = [theta y];
save Fig_1_gp0.d A -ascii -tabs
hold on
y = real(SSSTgradc (x, 1));
k = find (y ~= 0);
y(k) = -y(k)/min(y);
plot (theta, y)
A = [theta y];
save Fig_1_gp1.d A -ascii -tabs
y = real(SSSTgradc (x, 5));
k = find (y ~= 0);
y(k) = -y(k)/min(y);
plot (theta, y)
A = [theta y];
save Fig_1_gp5.d A -ascii -tabs
y = real(SSSTgradc (x, 25));
k = find (y ~= 0);
y(k) = -y(k)/min(y);
plot (theta, y)
A = [theta y];
save Fig_1_gp25.d A -ascii -tabs
y = real(SSSTgradc (x, 100));
k = find (y ~= 0);
y(k) = -y(k)/min(y);
plot (theta, y)
A = [theta -y/min(y)];
save Fig_1_gp100.d A -ascii -tabs
