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

% ============================================
% fill3dholes.m
%
% fill 3D holes
% connectivity 1=(6+,18), 2=(18,6+), 3=(6,26), 4=(26,6)
% epsilon: hole size, holes bigger than epsilon won't be filled
% 
% Li Shen
% 11/09/2004 - create

function obj = fill3dholes(obj,conn,epsilon,name)

d = size(obj); d = d+4;

% original binary image X
X = zeros(d);
X(3:end-2,3:end-2,3:end-2) = obj;
Xc = 1-X; % X compliment

% set distance map phi and hole size epsilon
phi = bwdist(X);
idx = find(phi>epsilon);
mask = ones(d); mask(idx) = 0;

% find bounding box B
idx = find(X~=0);
[xs,ys,zs] = ind2sub(d,idx);
xs = min(xs):max(xs);
ys = min(ys):max(ys);
zs = min(zs):max(zs);
B = zeros(d); B(xs,ys,zs) = 1;
B = B.*mask; % apply mask
Bc = 1-B; % B compliment

% initialize topological hull Y
Y = B;

% neighbour template
switch conn
    case {2,4} % (18,6+), (26,6)
		nbt(:,:,1) = [0 0 0; 0 1 0; 0 0 0];
		nbt(:,:,2) = [0 1 0; 1 0 1; 0 1 0];
		nbt(:,:,3) = [0 0 0; 0 1 0; 0 0 0];
    case {1} % (6+,18)
		nbt(:,:,1) = [0 1 0; 1 1 1; 0 1 0];
		nbt(:,:,2) = [1 1 1; 1 0 1; 1 1 1];
		nbt(:,:,3) = [0 1 0; 1 1 1; 0 1 0];
    case {3} % (6,26)
		nbt(:,:,1) = [0 1 0; 1 1 1; 0 1 0];
		nbt(:,:,2) = [1 1 1; 1 0 1; 1 1 1];
		nbt(:,:,3) = [0 1 0; 1 1 1; 0 1 0];
end

% initialize L list
L = [];
idxset = find(B-X==1); % B minus X
for i=1:length(idxset)
% 	disp(sprintf('%d/%d',i,length(idxset)));
    idx = idxset(i);
    [x,y,z] = ind2sub(d,idx);
%     nb = zeros(d); nb(x-1:x+1,y-1:y+1,z-1:z+1) = nbt;
%     tv = nb.*Bc; 
    tv = nbt.*Bc(x-1:x+1,y-1:y+1,z-1:z+1); 
    if sum(tv(:))>0
        L = [L,idx];
    end
end

dn = [3 3 3]; nb = ones(dn); nb(2,2,2) = 0;

% while loop
while ~isempty(L)
	disp(sprintf('|L|=%d',length(L)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	[sorted,order] = sort(phi(L));
% 	idx = L(order(end)); L = setdiff(L,idx);
    [xx,ord] = max(phi(L)); idx = L(ord); L(ord) = L(end); L = L(1:end-1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	[x,y,z] = ind2sub(d,idx);
	nbh = Y(x-1:x+1,y-1:y+1,z-1:z+1);
	issimple = simple([nbh(:)' conn]);
	if issimple==0
      % cannot remove the voxel
      continue;
	end
	% remove the voxel
	Y(idx) = 0; % remove idx from Y
	idxset = find(nb.*Xc(x-1:x+1,y-1:y+1,z-1:z+1).*Y(x-1:x+1,y-1:y+1,z-1:z+1)==1);
    if ~isempty(idxset)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 	    nb = zeros(d); nb(x-1:x+1,y-1:y+1,z-1:z+1) = 1; nb(x,y,z) = 0;
% 	    idxset = find(nb.*Xc.*Y==1);
        [xi,yi,zi] = ind2sub(dn,idxset);
        idxset = sub2ind(d,x-2+xi,y-2+yi,z-2+zi);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	    L = union(L,idxset);
    end
end

disp(sprintf('original %d, after_filling %d, difference %d',sum(X(:)),sum(Y(:)),length(find((Y-X)))));

figure;

subplot(2,2,1); [vertices, faces] = gensfdata(X);
patches = patch('faces', faces, 'vertices', vertices, ...
			'FaceVertexCData', ones(size(vertices,1),1)*[0.8 0.8 0.8], ...
			'FaceColor', 'interp', 'EdgeColor', 'none');
l1 = light('position', [-1 -1 .1], 'color', [.3 .3 .3]);
l2 = light('position', [1 -1 .1], 'color', [.3 .3 .3]);
l3 = light('position', [-.5 1 .1], 'color', [.3 .3 .3]);
l4 = light;
material([.3 .4 .2 10]); lighting phong;
axis equal; view(3); ax = axis; title(strrep(name,'_','-'));

subplot(2,2,2); [vertices, faces] = gensfdata(Y);
patches = patch('faces', faces, 'vertices', vertices, ...
			'FaceVertexCData', ones(size(vertices,1),1)*[0.8 0.8 0.8], ...
			'FaceColor', 'interp', 'EdgeColor', 'none');
l1 = light('position', [-1 -1 .1], 'color', [.3 .3 .3]);
l2 = light('position', [1 -1 .1], 'color', [.3 .3 .3]);
l3 = light('position', [-.5 1 .1], 'color', [.3 .3 .3]);
l4 = light;
material([.3 .4 .2 10]); lighting phong;
axis equal; view(3); axis(ax);

if ~isempty(find(Y-X==1))
	subplot(2,2,3); [vertices, faces] = gensfdata(Y-X);
	patches = patch('faces', faces, 'vertices', vertices, ...
				'FaceVertexCData', ones(size(vertices,1),1)*[0.8 0.8 0.8], ...
				'FaceColor', 'interp', 'EdgeColor', 'none');
	l1 = light('position', [-1 -1 .1], 'color', [.3 .3 .3]);
	l2 = light('position', [1 -1 .1], 'color', [.3 .3 .3]);
	l3 = light('position', [-.5 1 .1], 'color', [.3 .3 .3]);
	l4 = light;
	material([.3 .4 .2 10]); lighting phong;
	axis equal; view(3); axis(ax); 
end

obj = Y(3:end-2,3:end-2,3:end-2); 

return;

%
% generate surface data structure
%

function [vertices, faces] = gensfdata(bim)

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

disp(sprintf('%d = %d vertices - %d faces',size(vertices,1)-size(faces,1),size(vertices,1),size(faces,1)));

return;