function y = WB_lookup (x, p, c)
% WB_lookup  Spherical Surface Spline in Tension lookup (C mex version)
%
%   y = WB_lookup (x, p, c)
%
%   x is cos(theta)
%   p is tension
%   c is the solution for all x
%
% Returns the Green function for a spherical surface spline
% in tension, following Wessel and Becker [2008], via
% precalculated look-up tables.  Initialize by doing
%	c = SSSTc (x, p);
% for the full range of x of interest (finely sampled).
% If p == 0 then we return Parker's [1994] minimum
% curvature solution instead (dilog(sin^2(theta/2)))
