function [m, a] = polycenter(x,y)
% POLYCENTER  Compute center of mass and area of polygon
%
% [M, A] = POLYCENTER(X, Y)
%
%   M is a 2-vector with the coordinates of the center of mass of a
%   polygon.
%
%   A is the polygon area.
%
%   X, Y are vectors with the coordinates of the polygon vertices.

% Author: Ramon Casero
% Copyright © 2010 University of Oxford
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
error( nargchk( 2, 2, nargin, 'struct' ) );
error( nargoutchk( 0, 2, nargout, 'struct' ) );

% compute polygon area
a = polyarea( x, y );

% close polygon, if open
if (x(1) ~= x(end)) && (y(1) ~= y(end))
    x = x( [1:end 1] );
    y = y( [1:end 1] );
end

% compute x-coordinate of the centroid
m(1) = sum( ...
    ( x( 1:end-1 ) + x( 2:end ) ) .* ( ...
    x( 1:end-1 ) .* y( 2:end ) - x( 2:end ) .* y( 1:end-1 ) ) ...
    ) / 6 / a;

% compute y-coordinate of the centroid
m(2) = sum( ...
    ( y( 1:end-1 ) + y( 2:end ) ) .* ( ...
    x( 1:end-1 ) .* y( 2:end ) - x( 2:end ) .* y( 1:end-1 ) ) ...
    ) / 6 / a;
