function d = gdmatrix(lat, lon)
% GDMATRIX  Matrix of geodesic distances between points on a sphere, in
% degrees of arc
%
% D = gdmatrix(LAT, LON)
%
%   LAT, LON are vectors with the latitude and longitude coordinates of a
%   set of points on the surface of a sphere.
%
%   D is a square matrix where D(i,j) is the geodesic distance between the
%   i-th and j-th points.

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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

% check arguments
error(nargchk(2, 2, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% number of points
N = length(lat);
if (length(lon) ~= N)
    error('LAT and LON must have the same size')
end

% compute distance matrix row by row
d = zeros(N);
for I = 1:N
    d(I, :) = distance(lat(I), lon(I), lat, lon);
end
