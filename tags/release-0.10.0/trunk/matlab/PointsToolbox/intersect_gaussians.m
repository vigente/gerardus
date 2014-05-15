function x = intersect_gaussians(mu1, mu2, sigma1, sigma2)
% INTERSECT_GAUSSIANS  Intersection of two Gaussians (analytic solution)
%
% X = intersect_gaussians(MU1, MU2, SIGMA1, SIGMA2)
%
%   This function finds the intersection points of two Gaussian
%   distributions, with means MU1, MU2, and standard deviations SIGMA1,
%   SIGMA2, respectively. The expression for a Gaussian distribution is:
%
%     f(x) = 1/sigma*sqrt(2*pi) exp(-(x-mu)^2 / (2*sigma^2))
%
%   The solution is analytical, and it was computed using Mathematica (see
%   image intersect_gaussians.png in the current directory for the
%   solution).
%
%   MU1, MU2, SIGMA1, SIGMA2 can be scalars or column vectors of the same
%   length.
%
%   X is a two-column matrix, where each 

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2012 University of Oxford
% Version: 0.1.0
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
narginchk(4, 4);
nargoutchk(0, 1);

if ((size(mu1, 2) > 1) || (size(mu2, 2) > 1) ...
        || (size(sigma1, 2) > 1) || (size(sigma2, 2) > 1))
    error('MU1, MU2, SIGMA1, SIGMA2 must be scalars or column vectors')
end
if ((size(mu1, 1) ~= size(mu2, 1)) ...
        || (size(mu1, 1) ~= size(sigma1, 1)) ...
        || (size(mu1, 1) ~= size(sigma2, 1)))
    error('MU1, MU2, SIGMA1, SIGMA2 must have the same number of elements')
end

% init output
x = nan(size(mu1, 1), 2);

% identify identical gaussians, these return NaN, because all their points
% are identical
idx = (sigma1 == sigma2) & (mu1 == mu2);
if ~isempty(find(idx, 1))
    x(idx, :) = nan;
end

% identify pairs of gaussians with the same standard deviations
idx = (sigma1 == sigma2) & (mu1 ~= mu2);
if ~isempty(find(idx, 1))
    x(idx, 1) = (mu1(idx) + mu2(idx)) * .5;
    x(idx, 2) = Inf * (mu2(idx) - mu1(idx));
end

% no need to compute anything else if we only have equal sigmas
idx = (sigma1 ~= sigma2);
if ~isempty(find(idx, 1))
    
    % compute solutions when the sigmas are different
    mu1 = mu1(idx);
    mu2 = mu2(idx);
    sigma1 = sigma1(idx);
    sigma2 = sigma2(idx);
    SIGMA1 = sigma1 .* sigma1;
    SIGMA2 = sigma2 .* sigma2;
    
    aux1 = mu2 .* SIGMA1 - mu1 .* SIGMA2;
    aux2 = sigma1 .* sigma2 .* sqrt((...
        (mu1 - mu2).^2 + 2 * (SIGMA2 - SIGMA1) .* log(sigma2 ./ sigma1) ...
        ));
    aux3 = 1./(SIGMA1 - SIGMA2);
    
    % compute intersection points
    x(idx, 1) = (aux1 - aux2) .* aux3;
    x(idx, 2) = (aux1 + aux2) .* aux3;

end

% sort solutions from smaller to larger
x = sort(x, 2, 'ascend');
