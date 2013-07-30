function con2 = congrow(con)
% CONGROW  Grow the local neighbourhoods in a connectivity matrix
%
% CON2 = congrow(CON)
%
%   CON is a square connectivity matrix, typically sparse. CON(i,j)>0 means
%   that vertices i and j are connected with an edge.
%
%   CON2 is a connectivity matrix of the same size as CON. The connectivity
%   has been increased, by connecting each vertex to its neighbours'
%   neighbours. This makes the local neighbourhoods grow.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
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
narginchk(1, 1);
nargoutchk(0, 1);

if (size(con, 1) ~= size(con, 2))
    error('CON must be a square matrix')
end

% duplicate the connectivity matrix to add the new connections (it is
% important that we don't use the original matrix to store new connections,
% or else we can get running conditions where local neighbourhoods grow by
% more than one jump)
con2 = con;

% loop each vertex
for I = 1:size(con, 1)
    % neighbours of current vertex
    nei = con(:, I) > 0;
    
    % neighbours of the neighbours
    nei = sum(con(:, nei), 2) > 0;
    
    % connect current vertex to all its neighbours
    con2(I, nei) = 1;
    con2(nei, I) = 1;
end

% by convention, vertices are not connected to themselves
con2(1:size(con2, 1)+1:numel(con2)) = 0;
