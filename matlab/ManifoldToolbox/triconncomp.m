function [N, trilab, xlab] = triconncomp(tri, x)
% TRICONNCOMP  Find connected components in mesh
%
% [N, TRILAB, XLAB] = triconncomp(TRI, X)
%
%   Return the connected components found in the mesh described by TRI, X.
%
%   TRI is a 3-column matrix. Each row represents the indices of the three
%   vertices that form a triangle. TRI as a whole represents the closed
%   surface.
%
%   X is a 3-column matrix. Each row represents the Cartesian coordinates
%   of a vertex on the surface, indexed by TRI values.
%
%   N is the number of connected components found by the algorithm.
%
%   TRILAB and XLAB are vectors that label TRI and X, respectively. The
%   label indicates which connected components they belong to.

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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

% check arguments
narginchk(2, 2);
nargoutchk(0, 3);

% save triangle matrix for the labels
if (nargout > 1)
    tri0 = tri;
end

% sort triangles, so that the same vertices are nearby, and we don't need
% to create first lots of many separated connected components, and then
% merge them
tri = sort(tri, 2, 'ascend');
tri = sortrows(tri);

% init output
cc = {unique(tri(1, :))};

% assign vertices to separate connected components
for  I = 2:size(tri, 1)
    
    % find components that are incident to the current triangle
    idx = find(cellfun(@(x) any(intersect(tri(I, :), x)), cc));
    
    % if the triangle is in no component, then we have to create a new one
    if isempty(idx)
        cc{end+1} = tri(I, :);
        
        % if the triangle is incident to one component, we add the rest of
        % the vertices to it
    elseif (length(idx) == 1)
        cc{idx} = union(cc{idx}, tri(I, :));
        
        % if the triangle is incident to two components, we add the
        % triangle to one of them, and then combine them both
    elseif (length(idx) == 2)
        cc{idx(1)} = union(cc{idx(1)}, tri(I, :));
        cc{idx(1)} = union(cc{idx(1)}, cc{idx(2)});
        cc(idx(2)) = [];

        % if the triangle is incident to three components, we add the
        % triangle to one of them, and then combine all three of them
    elseif (length(idx) == 3)
        cc{idx(1)} = union(cc{idx(1)}, tri(I, :));
        cc{idx(1)} = union(cc{idx(1)}, cc{idx(2)});
        cc{idx(1)} = union(cc{idx(1)}, cc{idx(3)});
        cc(idx(2:3)) = [];
        
    else
        error('Assert fail: A triangle cannot have more than three vertices')
    end
    
end

% number of connected components
N = length(cc);

% if the user didn't ask for the labels, we don't need to waste time
% generating them
if (nargout < 3)
    return
end

% create labels for vertices
xlab = zeros(size(x, 1), 1);
for I = 1:N
    xlab(cc{I}) = I;
end

% recover the original ordering of the triangles
tri = tri0;

% create labels for triangles (the triangle will have the same label as its
% first vertex, as by construction all vertices in the triangle have the
% same label)
trilab = xlab(tri(:, 1));
