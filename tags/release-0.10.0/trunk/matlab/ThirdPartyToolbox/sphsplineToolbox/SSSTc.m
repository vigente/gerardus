function y = SSSTc (x, p)
% SSSTc  Spherical Surface Spline in Tension (C mex version)
%
%   y = SSSTc (x, p)
%
%   x is cos(theta) in range -1 <= x <= 1
%   p is tension (p >= 0)
%
% Returns the Green's' function for a spherical surface spline
% in tension, following Wessel and Becker [2008].
% If p == 0 or not given then we return Parker's [1994] minimum
% curvature solution instead.

% Nothing here - we execute the mex code
