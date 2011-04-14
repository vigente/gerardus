function [rforms, rformh, rformt] = pts_simil_map_params(y, x)
% PTS_SIMIL_MAP_PARAMS  Compute similarity transformation parameters
% between sets of points with unknown correspondence (Procrustes is used)
%
% [RFORMS, RFORMH, RFORMT] = PTS_SIMIL_MAP_PARAMS(Y, X)
%
%   In their simplest form:
%
%     Y: target point set
%     X: source point set
%
%   A point set is a (P,2)-matrix, where each row has the Cartesian
%   coordinates of a 2D point.
%
%   RFORMS, RFORMH and RFORMT are matrices with the parameters of the
%   similarity transformation that better maps X onto Y, in the least
%   squares sense.
%
%     RFORMS: s or scaling factor, (N,1)-vector
%     RFORMH: h or rotation matrix, (2,2,N)-matrix
%     RFORMT: t or translation vector, (N,2)-vector
%
%   The similarity transformation on a point [x1 x2] is
%
%     [z1 z2] = s * [x1 x2] * h + t
%
%   Input argument allowed sizes:
%
%     Y is (P,2), X is (P,2,N) => RFORMS is (N,1)-vector
%        RFORS(i) corresponds to X(:,:,i) -> Y
%        (similarly with RFORMH, RFORMT)
%
%     Y is (P,2,N), X is (P,2) => RFORMS is (N,1)-vector
%        RFORS(i) corresponds to X -> Y(:,:,i)
%        (similarly with RFORMH, RFORMT)
%
%     Y is (P,2,N), X is (P,2,N) => RFORMS is (N,1)-vector
%        RFORS(i) corresponds to X(:,:,i) -> Y(:,:,i)
%        (similarly with RFORMH, RFORMT)
%
% NB. The reason for not using a struct to hold all parameters is because
% it is e.g. ~20% faster to do the rotations with normal matrices than with
% matrices in structs.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2007-2011 University of Oxford
% Version: 1.0.0
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
error(nargoutchk(0, 3, nargout, 'struct'));

% get sizes
P = size(y, 1); % P: number of points in the configuration
if (P ~= size(x, 1))
    error('X and Y must have the same number of rows')
end
if (size(x, 2) ~= 2 || size(y, 2) ~= 2)
    error('Points in X or Y are not 2-dimensional')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Y is (P,2), X is (P,2,N) => RFORMS is (N,1)-vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (size(y, 3) == 1)
    N = size(x, 3); % N: number of source configurations
    
    % init output
    [rforms, rformh, rformt] = init_output(N);
    
    for I = 1:N
        
        % compute parameters
        [rforms(I), rformh(:, :, I), rformt(I, :)] = ...
            map_params(y, x(:, :, I));

    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Y is (P,2,N), X is (P,2) => RFORMS is (N,1)-vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (size(x, 3) == 1)
    N = size(y, 3); % N: number of target configurations
    
    % init output
    [rforms, rformh, rformt] = init_output(N);
    
    for I = 1:N

        % compute parameters
        [rforms(I), rformh(:, :, I), rformt(I, :)] = ...
            map_params(y(:, :, I), x);

    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Y is (P,2,N), X is (P,2,N) => RFORMS is (N,1)-vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    N = size(x, 3); % N: number of target configurations
    if (N ~= size(y, 3))
        error('Different number of point configurations in X and Y')
    end
    
    % init output
    [rforms, rformh, rformt] = init_output(N);
    
    for I = 1:N

        % compute parameters
        [rforms(I), rformh(:, :, I), rformt(I, :)] = ...
            map_params(y(:, :, I), x(:, :, I));
        
    end

end

% initialize output
function [rforms, rformh, rformt] = init_output(N)
rforms = zeros(N, 1);
rformh = zeros(2, 2, N);
rformt = zeros(N, 2);

% compute output
function [rforms, rformh, rformt] = map_params(y, x)
% map X onto Y
[~, ~, ~, s1, s2, rformh, t1, t2] = ...
    pts_procrustes(y, x);

% total scaling parameter
rforms = s2 / s1;

% compute composite rotation matrix
% already computed above

% total translation parameter
rformt = (t2 - t1) ./ s1;
