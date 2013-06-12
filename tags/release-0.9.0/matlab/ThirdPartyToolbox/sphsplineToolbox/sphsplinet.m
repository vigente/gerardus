function answer = sphsplinet(loni, lati, zi, lono, lato, p)
% SPHSPLINET  Gridding using a spherical surface spline in tension
%
% zo = sphsplinet (loni, lati, zi, lono, lato[, p])
%
% Here, loni, lati, zi are arrays with the input data constraints
% lono, lato are the desired output locations where the surface
% should be evaluated.
% p, if present, is the tension [Default is no tension]
% lons and lats are assumed to be in degrees

% $Id: sphsplinet.m,v 1.1.1.1 2008/05/09 21:34:52 myself Exp $
% P. Wessel, SOEST, U of Hawaii, April 2008 (pwessel@hawaii.edu)

if nargin == 5
    p = 0;
end

% Get Cartesian vectors
[x y z] = sph2cart (deg2rad(loni(:)), deg2rad(lati(:)), 1);
[xo yo zo] = sph2cart (deg2rad(lono(:)), deg2rad(lato(:)), 1);

% Find and remove mean z

zmean = mean (zi);
zrange = max(zi(:)) - min(zi(:));
zi = (zi - zmean) / zrange; % Normalize z data

% Set up linear system

n = length (zi);
A = zeros (n);
for row = 1:n
    g = ([x y z] * [x(row); y(row); z(row)])';
    A(row,:) = SSST(g,p);
end

% Solve for coefficients

f = A \ zi;

% Evaluate output solution

answer = zeros (size (lono));

for i = 1:length (xo)
    g = ([x y z] * [xo(i); yo(i); zo(i)])';
    answer(i) = zmean + zrange * SSST(g,p) * f;
end
