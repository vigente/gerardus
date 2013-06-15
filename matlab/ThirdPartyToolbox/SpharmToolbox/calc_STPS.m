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

function [new_ldmk2, new_verts] = calc_STPS(ldmk1, ldmk2, verts)

% Calculate c and d for calculating displacement in longitudinal and
% colatitudinal directions
[c_longitude, d_longitude, c_colatitude, d_colatitude] = calculate_displacement(ldmk1, ldmk2);

% Move spherical vertices to new locations
new_verts = move_spherical_vertices(verts, ldmk2, c_longitude, d_longitude, c_colatitude, d_colatitude);

% Calculate new landmark positions
new_ldmk2 = move_spherical_vertices(ldmk2, ldmk2, c_longitude, d_longitude, c_colatitude, d_colatitude);

return;

%
% move spherical vertices
%

function new_verts = move_spherical_vertices(vs, ldmks, c_longitude, d_longitude, c_colatitude, d_colatitude)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate new positions of sample points (vertices)
% Convert Cartesian coordinate into spherical coordinate
[longitude_vs, latitude_vs] = cart2sph(vs(:,1), vs(:,2), vs(:,3));

ind = find(longitude_vs<0); longitude_vs(ind) = longitude_vs(ind)+2*pi;
colatitude_vs = pi/2-latitude_vs;

R_vs = calculate_R_matrix(vs, ldmks);
R_vs = R_vs';

% Calculate displacement in longitude and colatitude directions
u_longitude_vs = c_longitude'*R_vs+d_longitude;
u_colatitude_vs = c_colatitude'*R_vs+d_colatitude;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate new spherical coordinates
new_longitude_vs = longitude_vs+u_longitude_vs';
new_colatitude_vs = colatitude_vs+u_colatitude_vs';

new_latitude_vs = pi/2-new_colatitude_vs;

% Limit the value of latitude between pi/2 and -pi/2
lg_LT=find(new_latitude_vs > (pi/2)); new_latitude_vs(lg_LT) = pi/2;
sm_LT=find(new_latitude_vs < (-1*pi/2)); new_latitude_vs(sm_LT) = -1*pi/2;

% Convert spherical coordinate into Cartesian coordinate
[x, y, z] = sph2cart(new_longitude_vs, new_latitude_vs, ones(size(new_longitude_vs,1),1));

new_verts = [x y z];

return;


%
% Calculate c vector and a scalar d in displacement mapping function u
%

function [c_longitude, d_longitude, c_colatitude, d_colatitude] = calculate_displacement(ldmk1, ldmk2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate spherical coordinates of two sets of landmarks in [longitude,
% colatitude]
% Spherical coordinates of landmark 1
[longitude1, latitude1] = cart2sph(ldmk1(:,1), ldmk1(:,2), ldmk1(:,3)); 
indp1 = find(longitude1<0); longitude1(indp1) = longitude1(indp1)+2*pi;
colatitude1 = pi/2-latitude1;

% Spherical coordinates of landmark 2
[longitude2, latitude2] = cart2sph(ldmk2(:,1), ldmk2(:,2), ldmk2(:,3));  
indp2 = find(longitude2<0); longitude2(indp2) = longitude2(indp2)+2*pi;
colatitude2 = pi/2-latitude2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate mapping function u for longitude and colatitude
% u_longitude(p) = displacement in longitude
% u_colatitude(p) = displacement in colatitude
% u = sum(c*R)+d

% First, calculate c and d from the known points (landmarks)
Rn = calculate_R_matrix(ldmk2,ldmk2);
iRn =  inv(Rn);

T = ones(size(ldmk2,1),1);
I = diag(ones(size(ldmk2,1), 1));

z_longitude = longitude1-longitude2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPORTANT : Minimize longituidinal displacement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lg_IDX = find(z_longitude >= pi);
z_longitude(lg_IDX) = z_longitude(lg_IDX)-2*pi;

sm_IDX = find(z_longitude <= -1*pi);
z_longitude(sm_IDX) = z_longitude(sm_IDX)+2*pi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z_colatitude = colatitude1-colatitude2;

d_longitude = (inv(T'*iRn*T)*T'*iRn)*z_longitude;
c_longitude = iRn*(z_longitude-d_longitude*T);

d_colatitude = (inv(T'*iRn*T)*T'*iRn)*z_colatitude;
c_colatitude = iRn*(z_colatitude-d_colatitude*T);

return;

%
% calculate R matrix
%

function R = calculate_R_matrix(p1, p2)
[x1, y1] = size(p1);
[x2, y2] = size(p2);

R = zeros(x1, x2);

Z = p1*(p2');  % Inner product matrix (cos(r(p1,p2)))
W = (1-Z)/2; 

ind = find(W==0);

Q = (log(1+1./sqrt(W)).*(12*(W.^2)-4*(W))-12*(W.^(3/2))+6*(W)+1)/2;  
Q(ind) = 0.5;
R = real((Q-1/3)/(4*pi));  

return;
