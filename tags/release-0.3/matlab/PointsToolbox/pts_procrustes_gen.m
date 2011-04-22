function [y, x, rho, sse] = pts_procrustes_gen(x, tol, VERB)
% PTS_PROCRUSTES_GEN Generalized Least-Squares Fit Orthogonal Procrustes
% Analysis
%
% [YM, Y] = PTS_PROCRUSTES_GEN(X)
%
%    X is a (P,K, N)-volume of sets of points, where P is the number of
%    points with dimension K, and N the number of sets.
%
%    YM is the mean set or consensus matrix.
%
%    Y is a volume with scaled and zero centered versions of X [1].
%    Besides, sets in Y has been rotated to match YM using Least-Squares
%    Fit Orthogonal Procrustes analysis (a form of similarity registration)
%    as described in [1, pp. 44-45].
%
%    An error in [1] that caused the algorithm to oscillate and not
%    converge in some cases has been fixed.
%
% [YM, Y] = PTS_PROCRUSTES_GEN(X, TOL, VERB)
%
%    The iterative algorithm stops when the improvement in SSE is < TOL. By
%    default, TOL = 1e-4.
%
%    VERB is a flag to make the algorithm verbose and show the evolution of
%    the SSE. By default, VERB=false.
%
% [YM, Y, RHO, SSE] = PTS_PROCRUSTES_GEN(...)
%
%    RHO is a vector with N elements, that are the scale factors for each
%    set. Note: rho is _not_ the size factor between the consensus and each
%    set of points. It is an internal parameter of the algorithm [1] that
%    scales each configuration with respect to the others.
%
%    SSE is the residual sum of squares of YM.
%
%
%  [1] Rohlf, F. & Slice, D. Extensions of the Procrustes method for the
%  optimal superimposition of landmarks. Systematic Zoology, 1990, 39,
%  40-59.
%
%  [2] Cootes, T. Image Processing and Analysis, Ch 7: Model-Based Methods
%  in Analysis of Biomedical Images Oxford University Press, 2000 ,
%  223-248.

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
error(nargchk(1, 3, nargin, 'struct'));
error(nargoutchk(0, 4, nargout, 'struct'));

% defaults
if (nargin < 2 || isempty(tol))
    tol = 1e-4;
end
if (nargin < 3 || isempty(VERB))
    VERB = false;
end

% get sizes
[~, K, N] = size(x);

% init variables
h = zeros(K, K, N);

% center and normalize all configurations independently
x = pts_cn(x);

% chose as consensus the first configuration
y = x(:, :, 1);

% rotate each configuration to fit consensus
for I = 1:N
    h(:, :, I) = pts_rotmat(y, x(:, :, I));
    x(:, :, I) = x(:, :, I) * h(:, :, I);
end

% update consensus matrix
y = mean(y, 3);

% compute initial SSE
sse0 = N * (1 - trace(y * y'));
sse = Inf;

% set all scale factors to 1
rho = ones(N, 1);

iter = 0;
while (abs(sse0 - sse) > tol)
    
    % update SSE
    sse0 = sse;
    
    % rotate and scale each configuration to fit current consensus
    for I = 1:N % loop configurations
        
        % compute rotation matrix
        h = pts_rotmat(y, x(:, :, I));
        
        % rotate configurations to consensus
        x(:, :, I) = x(:, :, I) * (h .* rho(I));
        
    end
        
    % update consensus
    y = mean(x, 3);
        
    for I = 1:N % loop configurations

        % compute scale ratio
        rhoratio = sqrt(abs(trace(x(:, :, I) * y') / ...
            trace(x(:, :, I) * x(:, :, I)') / ...
            trace(y * y')));

        % update configurations
        x(:, :, I) = x(:, :, I) * rhoratio;
        
        % update scale factor with scale ratio
        rho(I) = rho(I) * rhoratio;

    end

    % update consensus
    y = mean(x, 3);

    % compute residual sum of squares
    sse = N * (1 - trace(y * y'));

    if (VERB)
        iter = iter + 1;
    
        disp(['Iter: ' num2str(iter) '. |sse0 - sse| > tol?: ' ...
            num2str(abs(sse0 - sse)) ' > ' num2str(tol)])
    end
    
end

% normalize configurations  NORM(YM)=1
s = norm(y(:));
x = x / s;
y = y / s;
sse = N * (1 - trace(y * y'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % check Gower1975 condition as referenced by Rohlf and Slice 1990.
% function y = check(x)
% y = 0;
% for I = 1:size(x, 3)
%     y = y + trace(x(:, :, I) * x(:, :, I)');
% end
