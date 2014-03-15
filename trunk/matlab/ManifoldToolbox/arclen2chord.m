function d = arclen2chord(d, R)
% ARCLEN2CHORD  Sphere arc length distances to Euclidean chord length
% distances.
%
% DCHORD = arclen2chord(DARCLEN, R)
%
%   DARCLEN is a scalar, vector, matrix or general array with geodesic
%   (great circle) distances. Units are metres.
%
%   R is the radius of the sphere. Units are metres. By default, R=1.
%
%   DCHORD has the same size as DARCLEN, with the corresponding Euclidean
%   chord length distances. Units are metres.
%
% See also: chord2arclen.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
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
narginchk(1, 2);
nargoutchk(0, 1);

% defaults
if (nargin < 2 || isempty(R))
    R = 1.0;
end

% convert arc length distances to chord length distances (both in metres)
d = 2*R * sin(d/(2*R));
