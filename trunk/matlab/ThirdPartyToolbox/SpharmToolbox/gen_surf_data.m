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

%
% generate voxel surface from a binary image
%

function [vertices, faces] = gen_surf_data(bim,origin,vxsize)
vertices = []; faces = [];

if isempty(bim)
    return;
end
roi = bim; DIM = size(roi);

% make a work area so that all border voxels belong to the background
d = DIM+2; w = zeros(d);
w(2:d(1)-1,2:d(2)-1,2:d(3)-1) = roi;

% identifiy significant vertices (each vertice connects 8 adjacent voxels)
[xs, ys, zs] = meshgrid(1:d(1)-1,1:d(2)-1,1:d(3)-1);
len = prod(size(xs));
xs = reshape(xs,len,1);
ys = reshape(ys,len,1);
zs = reshape(zs,len,1);
nbsum = zeros(d);
inds = sub2ind(d,xs,ys,zs);
nbsum(inds) = sum([w(inds),...
                   w(sub2ind(d, xs+1,ys,  zs  )),...
                   w(sub2ind(d, xs,  ys+1,zs  )),...
                   w(sub2ind(d, xs+1,ys+1,zs  )),...
                   w(sub2ind(d, xs,  ys,  zs+1)),...
                   w(sub2ind(d, xs+1,ys,  zs+1)),...
                   w(sub2ind(d, xs,  ys+1,zs+1)),...
                   w(sub2ind(d, xs+1,ys+1,zs+1)),...
                   ]');
vertinds = find(nbsum~=0 & nbsum~=8);
[xs1,ys1,zs1] = ind2sub(d,vertinds); % significant ones
vnum = length(vertinds);
vertices = [xs1-1,ys1-1,zs1-1]; % remove the margin of the work area

% identify square faces on the surface
vertinds = sub2ind(DIM+1,xs1,ys1,zs1); % values of x, y, z could be 0
verts(1:prod(DIM+1)) = NaN;
verts(vertinds) = 1:vnum;

faces = [];

% front
tempind = find(w(inds)==1 & w(sub2ind(d, xs+1,ys,zs))==0); idx = inds(tempind);
[x,y,z] = ind2sub(d,idx); v1s = verts(sub2ind(DIM+1,x,y,z));
v2s = verts(sub2ind(DIM+1,x,y-1,z));
v3s = verts(sub2ind(DIM+1,x,y-1,z-1));
v4s = verts(sub2ind(DIM+1,x,y,z-1));
faces(end+1:end+length(x),:) = [v1s; v2s; v3s; v4s]'; 

% right
tempind = find(w(inds)==1 & w(sub2ind(d, xs,ys+1,zs))==0); idx = inds(tempind);
[x,y,z] = ind2sub(d,idx); v1s = verts(sub2ind(DIM+1,x,y,z));
v2s = verts(sub2ind(DIM+1,x,y,z-1));
v3s = verts(sub2ind(DIM+1,x-1,y,z-1));
v4s = verts(sub2ind(DIM+1,x-1,y,z));
faces(end+1:end+length(x),:) = [v1s; v2s; v3s; v4s]'; 

% top
tempind = find(w(inds)==1 & w(sub2ind(d, xs,ys,zs+1))==0); idx = inds(tempind);
[x,y,z] = ind2sub(d,idx); v1s = verts(sub2ind(DIM+1,x,y,z));
v2s = verts(sub2ind(DIM+1,x-1,y,z));
v3s = verts(sub2ind(DIM+1,x-1,y-1,z));
v4s = verts(sub2ind(DIM+1,x,y-1,z));
faces(end+1:end+length(x),:) = [v1s; v2s; v3s; v4s]'; 

% rear
tempind = find(w(inds)==0 & w(sub2ind(d, xs+1,ys,zs))==1); idx = inds(tempind);
[x,y,z] = ind2sub(d,idx); v1s = verts(sub2ind(DIM+1,x,y,z));
v2s = verts(sub2ind(DIM+1,x,y,z-1));
v3s = verts(sub2ind(DIM+1,x,y-1,z-1));
v4s = verts(sub2ind(DIM+1,x,y-1,z));
faces(end+1:end+length(x),:) = [v1s; v2s; v3s; v4s]'; 

% left
tempind = find(w(inds)==0 & w(sub2ind(d, xs,ys+1,zs))==1); idx = inds(tempind);
[x,y,z] = ind2sub(d,idx); v1s = verts(sub2ind(DIM+1,x,y,z));
v2s = verts(sub2ind(DIM+1,x-1,y,z));
v3s = verts(sub2ind(DIM+1,x-1,y,z-1));
v4s = verts(sub2ind(DIM+1,x,y,z-1));
faces(end+1:end+length(x),:) = [v1s; v2s; v3s; v4s]'; 

% bottom
tempind = find(w(inds)==0 & w(sub2ind(d, xs,ys,zs+1))==1); idx = inds(tempind);
[x,y,z] = ind2sub(d,idx); v1s = verts(sub2ind(DIM+1,x,y,z));
v2s = verts(sub2ind(DIM+1,x,y-1,z));
v3s = verts(sub2ind(DIM+1,x-1,y-1,z));
v4s = verts(sub2ind(DIM+1,x-1,y,z));
faces(end+1:end+length(x),:) = [v1s; v2s; v3s; v4s]'; 

% scaling and translation
vertices = vertices.*vxsize(ones(1,size(vertices,1)),:);
vertices = vertices + origin(ones(1,size(vertices,1)),:);

return;