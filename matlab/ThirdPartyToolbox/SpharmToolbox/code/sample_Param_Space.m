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
% Goal: 
%   1. Create homogeneous distribution of sampling points on a unit sphere
%   2. Regular mesh sampling
%
% Li Shen 
% 10/14/2002 - create
% 11/03/2008 - modified by Sungeun Kim

function [vs, fs] = sample_Param_Space(res)

if res<=0  
    % homogeneour sampling using icosahedron
	[vs, fs] = icosahedron;
    vnum = size(vs,1); fnum = size(fs,1);
%    name = sprintf('subdiv_v%03d_f%04d',vnum,fnum);
    
	for i = 1:-res
        [vs, fs] = subdivision(vs, fs);
        vnum = size(vs,1); fnum = size(fs,1);
%        name = sprintf('subdiv_v%03d_f%04d',vnum,fnum);
	end
else
    [x, y, z] = sphere(res);
    [fs,vs,c] = surf2patch(x,y,z);

end

return;

%
% create icosahedron with 20 faces, 12 vertices, and 30 edges
%
function [vs, fs] = icosahedron

% create icosahedron
t1=3*pi/10;t=pi/5;r=.5*sec(t1);s=.5*tan(t1);a=.5*sqrt(.75-(r-s)*(r-s));b=sqrt(1-r*r);
ind2=[(1:2:11),(2:2:10)];ind3=[(11:-2:1),(10:-2:2)];t2=t*(0:2:20)';t3=t*(1:-2:-19)';
x=zeros(11,4);x(ind2,2)=r*cos(t2);x(ind3,3)=r*cos(t3);
y=zeros(11,4);y(ind2,2)=r*sin(t2);y(ind3,3)=r*sin(t3);
z=zeros(11,4);z(:,1)=(a+b)*ones(11,1);z(:,2)=a*ones(11,1);z(:,3)=-z(:,2);z(:,4)=-z(:,1);

vs(1,:)    = [x(1,1)    y(1,1)    z(1,1)];
vs(2:6,:)  = [x(1:2:10,2) y(1:2:10,2) z(1:2:10,2)];
vs(7:11,:) = [x(1:2:10,3) y(1:2:10,3) z(1:2:10,3)];
vs(12,:)   = [x(1,4)    y(1,4)    z(1,4)]; 

fs(1:5,:)  = [1 2 3; 1 3 4; 1 4 5; 1 5 6; 1 6 2];
fs(6:10,:) = [7 3 2; 8 4 3; 9 5 4; 10 6 5; 11 2 6];
fs(11:15,:) = [3 7 8; 4 8 9; 5 9 10; 6 10 11; 2 11 7];
fs(16:20,:) = [12 8 7; 12 9 8; 12 10 9; 12 11 10; 12 7 11];

len = sqrt(sum((vs.^2)'))';
vs = vs./len(:,[1 1 1]);

return;

%
% do subdivision
%
function [new_vs, new_fs] = subdivision(vs, fs)

vertnum = size(vs,1); facenum = size(fs,1);

% collect edges
edges = [fs(:,1) fs(:,2); fs(:,2) fs(:,3); fs(:,3) fs(:,1)];
ind = find(edges(:,1)<edges(:,2)); edges = edges(ind,:); edgenum = size(edges,1);

% generate new vertices
new_vs = (vs(edges(:,1),:)+vs(edges(:,2),:))/2;
len = sqrt(sum((new_vs.^2)'))';
new_vs = new_vs./len(:,[1 1 1]);

% generate faces and vertices
edges = edges(:,1) + j*edges(:,2);
for i=1:facenum
    es = min(fs(i,1), fs(i,2)) + j*max(fs(i,1), fs(i,2));
    vind_a(i,1) = find(edges==es) + vertnum;
    es = min(fs(i,2), fs(i,3)) + j*max(fs(i,2), fs(i,3));
    vind_b(i,1) = find(edges==es) + vertnum;
    es = min(fs(i,3), fs(i,1)) + j*max(fs(i,3), fs(i,1));
    vind_c(i,1) = find(edges==es) + vertnum;
end

new_fs = [fs(:,1) vind_a  vind_c;
          vind_a  fs(:,2) vind_b;
          vind_c  vind_b  fs(:,3);
          vind_a  vind_b  vind_c];
new_vs = [vs; new_vs];

return;