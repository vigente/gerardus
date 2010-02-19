function h = plotaxes(a, m, type, col)
% PLOTAXES  Plot axes of 3D coordinate system
%
% PLOTAXES(A, M, TYPE, COL)
%
%   A is a (3,3)-matrix with 3 column vectors. Each vector gives the
%   direction of one axis of a Coordinate system.
%
%   M is a 3-vector with the origin of coordinates. By default, M=[0,0,0].
%
%   TYPE is a cell array with the line types to plot (e.g. '--', ':'... see
%   "help plot" for more details). By default, type='-', i.e. a solid line.
%
%   COL is a cell array with the colours of each axis. By default,
%   col={'r', 'b', 'g'}. If only one string is provided, all axes are
%   displayed with the same colour.

% Copyright © 2010 University of Oxford
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
error( nargchk( 1, 4, nargin, 'struct' ) );
error( nargoutchk( 0, 1, nargout, 'struct' ) );

% defaults
if ( nargin < 2 || isempty( m ) )
    m = [0 0 0];
end
if ( nargin < 3 || isempty( type ) )
    type = {'-', '-', '-'};
end
if ( nargin < 4 || isempty( col ) )
    col = {'r', 'b', 'g'};
end
if ( ~iscell( type ) || length( type ) == 1 )
    if iscell( type )
        type = type{:};
    end
    type = { type, type, type };
end
if ( ~iscell( col ) || length( col ) == 1 )
    if iscell( col )
        col = col{:};
    end
    col = { col, col, col };
end

% save hold state for later
PHOLD = ishold;

% plot axes
for I = 1:size( a, 2 )
    h = plot3( [m(1), a(1, I)], [m(2), a(2, I)], [m(3), a(3, I)], [ type{I} col{I} ] );
    hold on
end

% recovert hold state
if PHOLD
    hold on
else
    hold off
end
