function [total_area, areas, sides] = trifacet_area3D(faces, vertices)
% trifacet_area3D 
%
% Function to calculate the surface area of a 3D mesh of triangles. This
% function uses Heron's forumla to calculate the sides of each triangle.
% 
% [TOTAL_AREA, AREAS, SIDES] = trifacet_area3D(FACES, VERTICES)
%  
%   FACES is a 3-column matrix. Each row of FACES contains the indices of
%   one triangle in the mesh.
%   VERTICES is a 3-column matrix (x,y and z).
%   Each row contains the x, y and z coordinates of one of the triangles.
%
%   TOTAL_AREA is the sum of all the triangle areas.
%   AREAS is a vector containing the areas of all the the
%   triangles.
%   SIDES is a 3-column matrix, containing the side lengths for each
%   triangle. The number of rows is the number of triangles.
%    

% Author: Jack Allen <jack.allen@jesus.ox.ac.uk>
% Adapted from: 
%    TRIFACET_SIGNED_AREA  (Gerardus)
%    Author: Ramon Casero <rcasero@gmail.com> 
%    Copyright © 2013 University of Oxford 
%    Version: 0.2.0 
% and
%   (side length and area calculations: Sean de Wolski
%   (http://www.mathworks.com/matlabcentral/answers/14928-area-of-triangle, accessed: May 2015))

% error messages
%(taken from 'trifacet_signed_area.m' (Author: Ramon Casero) in gerardus)
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