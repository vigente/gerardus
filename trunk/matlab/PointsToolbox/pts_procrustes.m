function [y1, y2, sse, s1, s2, h, t1, t2] = pts_procrustes(x1, x2)
% PTS_PROCRUSTES Least-Squares Fit Orthogonal Procrustes Analysis between
% two sets of points
%
% [Y1, Y2] = PTS_PROCRUSTES(X1, X2)
%
%    X1, X2 are sets of points with dimensions (P,K), where P is the
%    number of points in a set, and K the dimension of each point.
%
%    Y1, Y2 are scaled (square root of sum of squared distances = 1) and
%    zero centered versions of X1, X2. Besides, Y2 has been rotated to
%    match Y1 using Least-Squares Fit Orthogonal Procrustes analysis (a
%    form of similarity registration) as described in [1, pp. 41-42].
%
% [Y1, Y2, SSE, S1, S2, H, T1, T2] = PTS_PROCRUSTES(...)
%
%    SSE is the sum of squared differences between Y1 and Y2.
%
%    The similarity transformation that maps each point [x1, x2] to point 
%    [y1, y2] is formulated as
%
%       y1 = x1 .* S1 + T1
%       y2 = x2 * H .* S2 + T2
%
%    where 
%       H is the rotation matrix. Note: In 2D, the counterclockwise
%       definition of H is
%          H = [cos(phi) sin(phi) ; -sin(phi) cos(phi)]
%       and the clockwise
%          H = [cos(phi) -sin(phi) ; sin(phi) cos(phi)]
%       S1, S2 are the scaling factor, or sum of the squared distances of
%       each point to the object's centroid;
%       T1, T2 are translations expressed as (1, 2) vectors.
%
%  [1] Rohlf, F. & Slice, D. Extensions of the Procrustes method for the
%  optimal superimposition of landmarks. Systematic Zoology, 1990, 39,
%  40-59.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 1.0.0
% $Rev$
% $Date$
% 
% University of Oxford means the Chancellor, Masters and Scholars of
% the University of Oxford, having an administrative office at
% Wellington Square, Oxford OX1 2JD, UK. 
%
% This file is part of Gerardus.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. The offer of this
% program under the terms of the License is subject to the License
% being interpreted in accordance with English Law and subject to any
% action against the University of Oxford being under the jurisdiction
% of the English Courts.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% check arguments
error(nargchk(2, 2, nargin, 'struct'));
error(nargoutchk(0, 8, nargout, 'struct'));

% get sizes
[P, K] = size(x1);
if (P ~= size(x2, 1) || K ~= size(x2, 2))
    error('X1 and X2 must have the same dimensions')
end

% center and scale configurations
[y1, s1] = pts_cn(x1);
[y2, s2] = pts_cn(x2);

% compute rotation matrix
h = pts_rotmat(y1, y2);

% rotate second configuration
y2 = y2 * h;

% compute SSE of the match
sse = trace((y1 - y2) * (y1 - y2)');

% express translation operation in an explicit way
if (nargout > 6)
    t1 = y1(1, :) - x1(1, :) .* s1;
    t2 = y2(1, :) - (x2(1, :) * h) .* s2;
end
