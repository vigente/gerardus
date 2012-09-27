% test_intersect_gaussians.m

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

% x-axis values
x = linspace(-10, 10, 1001);

%% case where standard deviations are equal
mu1 = 0;
mu2 = -1;
sigma1 = 2;
sigma2 = 2;

% compute gaussians
f1 = 1/(sqrt(2*pi)*sigma1) * exp(-(((x-mu1)./sigma1).^2)/2);
f2 = 1/(sqrt(2*pi)*sigma2) * exp(-(((x-mu2)./sigma2).^2)/2);

% compute intersections
xint = intersect_gaussians(mu1, mu2, sigma1, sigma2);

% plot curve
hold off
plot(x, f1)
hold on
plot(x, f2, 'r')
for I = 1:2
    if ~isinf(xint(I))
        plot(xint([I I]), [0, max([f1(:); f2(:)])], 'g')
    end
end


%% case where standard deviations are different
mu1 = 0;
mu2 = 1;
sigma1 = 2;
sigma2 = 1;

% compute gaussians
f1 = 1/(sqrt(2*pi)*sigma1) * exp(-(((x-mu1)./sigma1).^2)/2);
f2 = 1/(sqrt(2*pi)*sigma2) * exp(-(((x-mu2)./sigma2).^2)/2);

% compute intersections
xint = intersect_gaussians(mu1, mu2, sigma1, sigma2);

% plot curve
hold off
plot(x, f1)
hold on
plot(x, f2, 'r')
for I = 1:2
    if ~isinf(xint(I))
        plot(xint([I I]), [0, max([f1(:); f2(:)])], 'g')
    end
end
