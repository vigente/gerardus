function w = pts_tps_weights(s, t)
% PTS_TPS_WEIGHTS  Compute weights and affine parameters of thin-plate
% spline warp for N-dimensional points
%
% W = PTS_TPS_WEIGHTS(S, T)
%
%    Compute the thin-plate splines warp given by [1], but with the
%    definition of the energy function in [2], both generalised from 2D to
%    N-dimensions.
%
%    For example, in 2D the equation to interpolate a point (x,y) given a
%    configuration of 2D points {Pi} is
%
%       f(x,y,z) = a1 + ax*x + ay*y + sum(wi*U(|Pi - (x,y)|))
%                                           i
%
%    W = [w1_x, ..., wN_x, a1_x, ax_x, ay_x ;
%          w1_y, ..., wN_y, a1_y, ax_y, ay_y]'
%
%    U(r) = r^2 log10(r)
%
%    S, T are matrices with the source and target point coordinates. They
%    do not need to have the same number of columns. 
%
%    For instance, if you want to interpolate the heights on a 2D map of
%    Earth, you would have 2 columns in S for the (x,y)-coordinates, and 1
%    column in T for the height values.
%
%    Another example, if you want to warp a 2D image, then you would have 2
%    columns both in S and T, with the coordinates of corresponding points.
%
%    You can also concatenate S matrices, and T matrices, and use the
%    resulting volumes as the input arguments. In that case, the weights
%    are computed independently for each pair of configurations, i.e.
%
%       W(:,:,1) for S(:,:,1) and T(:,:,1)
%       W(:,:,2) for S(:,:,2) and T(:,:,2)
%       etc.
%
% See also: pts_tps_map.
%
% [1] F.L. Bookstein. Principal Warps: Thin-plate splines and the
% decomposition of deformations. IEEE Trans on PAMI, 11(6), Jun 1989.
%
% [2] Barrodale et al. Warping digital images using thin plate splines.
% Pattern Recognition, 26(2):375--376, 1993.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2006-2013 University of Oxford
% Version: 0.5.0
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
narginchk(2, 2) ;
nargoutchk(0, 1) ;

% point dimensionality
D = size(s, 2);

% get sizes
P = size(s, 1); % P: number of points
if (P ~= size(t, 1))
    error('S and T must have the same number of points')
end
N = size(s, 3); % N: number of point configurations
if (N ~= size(t, 3))
    error('There must be the same number of point configurations in S and T')
end

% init output
w = zeros(P+D+1, size(t, 2), N);

% loop pairs of S and T point configurations
for I = 1:N % loop configurations

    % compute distances between source points without having to loop
    % between each pair of points in the configuration
    
    % for each dimension, first replicate the point coordinates, as an
    % intermediate result...
    %
    % For example, if s is a configuration of 5 points in 3D
    %
    %  s =
    %
    %  3.5784   -0.0631    1.4090
    %  2.7694    0.7147    1.4172
    % -1.3499   -0.2050    0.6715
    %  3.0349   -0.1241   -1.2075
    %  0.7254    1.4897    0.7172
    %
    % then
    %
    %  aux(:,:,1) =
    %
    %  3.5784    3.5784    3.5784    3.5784    3.5784
    %  2.7694    2.7694    2.7694    2.7694    2.7694
    % -1.3499   -1.3499   -1.3499   -1.3499   -1.3499
    %  3.0349    3.0349    3.0349    3.0349    3.0349
    %  0.7254    0.7254    0.7254    0.7254    0.7254
    %
    %
    %  aux(:,:,2) =
    %
    % -0.0631   -0.0631   -0.0631   -0.0631   -0.0631
    %  0.7147    0.7147    0.7147    0.7147    0.7147
    % -0.2050   -0.2050   -0.2050   -0.2050   -0.2050
    % -0.1241   -0.1241   -0.1241   -0.1241   -0.1241
    %  1.4897    1.4897    1.4897    1.4897    1.4897
    %
    % etc.

    aux = repmat(reshape(s(:, :, I), [P, 1, D]), ...
        [1, P, 1]);
    
    % ... and then substract the transpose of each matrix, square and sum
    %
    % K is a matrix where element (a,b) is the squared distance between
    % points a and b in the configuration
    %
    % the permute means: turn rows into columns, and viceversa, but leave
    % the 3rd component of aux (dimension) alone
    K = sum((aux - permute(aux, [2, 1, 3])) .^ 2, 3);

    % compute thin-plate energy function
    % U(r) = r^2 ln(r^2)
    warning('off', 'MATLAB:log:logOfZero');
    K = K .* log(K);
    K(isnan(K)) = 0;
    warning('on', 'MATLAB:log:logOfZero');

    % compute system matrix
    L = [K           ones(P, 1) s(:, :, I) ; ...
        ones(1, P)  zeros(1, D+1) ; ...
        s(:, :, I)' zeros(D, D+1)];

    % compute solution weights and affine parameters
    w(:, :, I) = L \ [t(:, :, I) ; zeros(D+1, size(t, 2))];

end
