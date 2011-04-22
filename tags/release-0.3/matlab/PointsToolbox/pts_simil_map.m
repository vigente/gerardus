function y = pts_simil_map(x, rforms, rformh, rformt)
% PTS_SIMIL_MAP  Apply similarity transformation to points
%
% Y = PTS_SIMIL_MAP(X, RFORMS, RFORMH, RFORMT)
%
%   X, Y are (P,2,N)-matrices of point configurations, i.e. each X(:,:,i)
%   point configuration is mapped to a Y(:,:,i) point configuration.
%
%   RFORMS, RFORMH and RFORMT are matrices with the parameters of the
%   similarity transformation that better maps X onto Y, in the least
%   squares sense.
%
%     RFORMS: s or scaling factor, (N,1)-vector
%     RFORMH: h or rotation matrix, (2,2,N)-volume
%     RFORMT: t or translation vector, (N,2)-matrix
%
%   The similarity transformation on a point [x1 x2] is
%
%     [y1 y2] = s * [x1 x2] * h + t
%
%   Input argument allowed sizes:
%
%     X is (P,2,N), RFORMS is scalar
%        (similarly with RFORMH, RFORMT)
%
%     X is (P,2), RFORMS is (N,1)-vector
%        (similarly with RFORMH, RFORMT)
%
%     X is (P,2,N), RFORMS is (N,1)-vector
%        (similarly with RFORMH, RFORMT)
%
% Y = PTS_SIMIL_MAP(X, Q)
%
%     Q is a 4-vector with the similarity transformation parameters
%
%     RFORMT = [Q(1), Q(2)]
%     RFORMS = 1+Q(3)
%     RFORMH = [cos(Q(4)) sin(Q(4));
%               -sin(Q(4)) cos(Q(4))]
%
%     Q can also be a matrix where each column is a 4-vector.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
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
error(nargchk(2, 4, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y = PTS_RIGID_MAP(X, Q)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin == 2)
    
    % reformat transformation parameters
    q = rforms;
    rformt = q(1:2, :)';
    rforms = q(3, :)' + 1;
    aux = reshape(q(4, :), 1, 1, size(q, 2));
    rformh = [cos(aux) sin(aux) ; ...
        -sin(aux) cos(aux)];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Y = PTS_RIGID_MAP(X, RFORMS, RFORMH, RFORMT)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get sizes
Nmap = size(rforms, 1); % Nmap: number of transformation
if (Nmap ~= size(rformh, 3) || Nmap ~= size(rformt, 1))
    error('RFORMS, RFORMH and RFORMT must have the same number of transformation parameters')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     X is (P,2,N), RFORMS is scalar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isscalar(rforms)
    
    if isstruct(x)
        
        % number of points in the point configuration
        Pin = size(x(1).in , 1);
        Pendo = size(x(1).endo , 1);
        Pepi = size(x(1).epi , 1);
        
        N = length(x); % N: number of point configurations
        
        % init output
        y(N) = struct(...
            'in', zeros(size(x(1).in)), ...
            'endo', zeros(size(x(1).endo)), ...
            'epi', zeros(size(x(1).epi)));

        % cache rformt
        rformt_in = ones(Pin, 1) * rformt;
        rformt_endo = ones(Pendo, 1) * rformt;
        rformt_epi = ones(Pepi, 1) * rformt;
        
        % apply rigid transformation to input
        for I = 1:N
            
            y(I).in = rforms * x(I).in * rformh + rformt_in;
            y(I).endo = rforms * x(I).endo * rformh + rformt_endo;
            y(I).epi = rforms * x(I).epi * rformh + rformt_epi;

        end

    else % ~isstruct(x)

        P = size(x, 1); % P: number of points in the point configuration
        N = size(x, 3); % N: number of point configurations

        % init output
        y = zeros(size(x));

        % cache rformt
        rformt = ones(P, 1) * rformt;

        % apply rigid transformation to input
        for I = 1:N

            y(:, :, I) = rforms * x(:, :, I) * rformh ...
                + rformt;

        end

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     X is (P,2), RFORMS is (N,1)-vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif (isstruct(x) && length(x) == 1)
    
    % number of points in the point configuration
    Pin = size(x(1).in , 1);
    Pendo = size(x(1).endo , 1);
    Pepi = size(x(1).epi , 1);

    % init output
    y(Nmap) = struct(...
        'in', zeros(size(x.in)), ...
        'endo', zeros(size(x.endo)), ...
        'epi', zeros(size(x.epi)));

    % apply rigid transformation to input
    for I = 1:Nmap

        y(I).in = rforms(I) * x.in * rformh(:, :, I) ...
            + ones(Pin, 1) * rformt(I, :);
        y(I).endo = rforms(I) * x.endo * rformh(:, :, I) ...
            + ones(Pendo, 1) * rformt(I, :);
        y(I).epi = rforms(I) * x.epi * rformh(:, :, I) ...
            + ones(Pepi, 1) * rformt(I, :);

    end
    
elseif (~isstruct(x) && size(x, 3) == 1)
    
    P = size(x, 1); % P: number of points in the point configuration

    % init output
    y = zeros([size(x) Nmap]);
    
    % apply rigid transformation to input
    for I = 1:Nmap

        y(:, :, I) = rforms(I) * x * rformh(:, :, I) ...
            + ones(P, 1) * rformt(I, :);

    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     X is (P,2,N), RFORMS is (N,1)-vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif isstruct(x)
    
    N = length(x); % N: number of point configurations
    
    % number of points in the point configuration
    Pin = size(x(1).in , 1);
    Pendo = size(x(1).endo , 1);
    Pepi = size(x(1).epi , 1);

    % init output
    y(N) = struct(...
        'in', zeros(size(x(1).in)), ...
        'endo', zeros(size(x(1).endo)), ...
        'epi', zeros(size(x(1).epi)));
    
    % apply rigid transformation to input
    for I = 1:N

        y(I).in = rforms(I) * x(I).in * rformh(:, :, I) ...
            + ones(Pin, 1) * rformt(I, :);
        y(I).endo = rforms(I) * x(I).endo * rformh(:, :, I) ...
            + ones(Pendo, 1) * rformt(I, :);
        y(I).epi = rforms(I) * x(I).epi * rformh(:, :, I) ...
            + ones(Pepi, 1) * rformt(I, :);

    end

else

    N = size(x, 3); % N: number of point configurations
    P = size(x, 1); % P: number of points in the point configuration
    
    % init output
    y = zeros(size(x));
    
    % apply rigid transformation to input
    for I = 1:N

        y(:, :, I) = rforms(I) * x(:, :, I) * rformh(:, :, I) ...
            + ones(P, 1) * rformt(I, :);

    end

end
