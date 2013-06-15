%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherical Harmonic Modeling and Analysis Toolkit (SPHARM-MAT) is a 3D 
% shape modeling and analysis toolkit. 
% It is a software package developed at Shenlab in Center for Neuroimaging, 
% Indiana University (SpharmMat@gmail.com, http://www.iupui.edu/~shenlab/)
% It is available to the scientific community as copyright freeware 
% under the terms of the GNU General Public Licence.
% 
% Copyright 2009, 2010, ShenLab, Center for Neuroimaging, Indiana University
% 
% This file is part of SPHARM-MAT.
% 
% SPHARM-MAT is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% SPHARM-MAT is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with SPHARM-MAT. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [areas, count] = calc_triangle_areas(verts, faces, space, focus)

if strcmp(space, 'object')
    tri_areas = cal_obj_area(verts, faces);
elseif strcmp(space, 'parameter')
    tri_areas = cal_par_area(verts, faces);
end

if strcmp(focus, 'vertex')
    areas = zeros(length(verts),1);
    count = zeros(length(verts),1);

% find all incident triangles upon each vertex and accumulate the areas of the triangles   
    for i=1:size(faces,2)
        [u1, m1, n1] = unique(faces(:,i));
        areas(u1) = areas(u1) + accumarray(n1, tri_areas);
        count(u1) = count(u1) + accumarray(n1, ones(length(n1),1));
    end
    
elseif strcmp(focus, 'triangle')
    areas = tri_areas;
    count = ones(length(tri_areas),1);
end

return;


%
% calculate relative areas of triangles on object surface net
%

function obj_area = cal_obj_area(vertices,faces)

A = faces(:,1); B = faces(:,2); C = faces(:,3);
a = sqrt(sum(((vertices(A,:)-vertices(B,:)).^(2))'))';
b = sqrt(sum(((vertices(B,:)-vertices(C,:)).^(2))'))';
c = sqrt(sum(((vertices(C,:)-vertices(A,:)).^(2))'))';
s = (a+b+c)/2;
obj_area = sqrt(s.*(s-a).*(s-b).*(s-c));
obj_area = obj_area/sum(obj_area);

return;


%
% calculate relative areas of spherical triangles in parameter space
%

function par_area = cal_par_area(vs,faces)

angles = [];
for j = 1:3
    % note that the order of A B C is clockwise (see 08-22-02.htm notes)
    A = vs(faces(:,j),:);
    B = vs(faces(:,mod(j,3)+1),:);
    C = vs(faces(:,mod(j-2,3)+1),:);
    y = A(:,1).*B(:,2).*C(:,3) - A(:,1).*B(:,3).*C(:,2) + ...
        A(:,2).*B(:,3).*C(:,1) - A(:,2).*B(:,1).*C(:,3) + ...
        A(:,3).*B(:,1).*C(:,2) - A(:,3).*B(:,2).*C(:,1);
    x = B(:,1).*C(:,1) + B(:,2).*C(:,2) + B(:,3).*C(:,3) - ...
       (A(:,1).*C(:,1) + A(:,2).*C(:,2) + A(:,3).*C(:,3)).* ...
       (A(:,1).*B(:,1) + A(:,2).*B(:,2) + A(:,3).*B(:,3));
    angles(:,j) = atan2(y,x); 
end
ind = find(angles<0);
angles(ind) = angles(ind) + 2*pi;
par_area = sum(angles')' - pi;
par_area = par_area/(4*pi);

return;
