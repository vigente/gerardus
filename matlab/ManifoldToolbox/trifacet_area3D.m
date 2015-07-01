function [total_area, areas, sides] = trifacet_area3D(faces, vertices)
% TRIFACET_AREA3D  Surface area of triangles in triangular mesh.
% 
% Function to calculate the surface area of a 3D mesh of triangles. This
% function uses Heron's forumla to calculate the sides of each triangle.
% 
% [TOTAL_AREA, AREAS, SIDES] = TRIFACET_AREA3D(FACES, VERTICES)
%  
%   FACES is a 3-column matrix. Each row of FACES contains the indices of
%   one triangle in the mesh.
%
%   VERTICES is a 3-column matrix (x,y and z).
%   Each row contains the x, y and z coordinates of one of the triangles.
%
%   TOTAL_AREA is the sum of all the triangle areas.
%
%   AREAS is a vector containing the areas of all the the
%   triangles.
%
%   SIDES is a 3-column matrix, containing the side lengths for each
%   triangle. The number of rows is the number of triangles.
%
% See also: trifacet_signed_area.

% Author: Jack Allen <jack.allen@jesus.ox.ac.uk>
% Copyright © 2015 University of Oxford
% Version: 0.1.1
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

if (size(faces, 2) ~= 3)
    error('TRI must have 3 columns')
end
if (size(vertices, 2) ~= 3)
    error('vertices must have 3 columns')
end

%allocate matrix to store side lengths of triangles
[nTrifacRows, ~] = size(faces);
sides = zeros(nTrifacRows,3);
%find side lengths for each triangle
%side a
sides(:,1) = sqrt(...
    (vertices(faces(:,1),1) - vertices(faces(:,2),1)).^2 + ... % sqrt((x1 - x2)^2) 
    (vertices(faces(:,1),2) - vertices(faces(:,2),2)).^2 + ... % sqrt((y1 - y2)^2)
    (vertices(faces(:,1),3) - vertices(faces(:,2),3)).^2); % sqrt((z1 - z2)^2)
%side b
sides(:,2) = sqrt(...
    (vertices(faces(:,2),1) - vertices(faces(:,3),1)).^2 + ...
    (vertices(faces(:,2),2) - vertices(faces(:,3),2)).^2 + ...
    (vertices(faces(:,2),3) - vertices(faces(:,3),3)).^2);
%side c
sides(:,3) = sqrt(...
    (vertices(faces(:,1),1) - vertices(faces(:,3),1)).^2 + ...
    (vertices(faces(:,1),2) - vertices(faces(:,3),2)).^2 + ...
    (vertices(faces(:,1),3) - vertices(faces(:,3),3)).^2);



% semi-perimeters s
s = (sides(:,1) + sides(:,2) + sides(:,3))./2;

areas = sqrt(s.*(s-sides(:,1)).*(s-sides(:,2)).*(s-sides(:,3)));

% total surface area of 3D mesh (scalar)
total_area = sum(sqrt(s.*(s-sides(:,1)).*(s-sides(:,2)).*(s-sides(:,3))));


end
