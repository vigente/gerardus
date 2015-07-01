function [intensities1, intensities2,intersectCoords] =  scimat_lineintersection(scimat1, scimat2, coords1, coords2, res)

% SCIMAT_LINEINTERSECTION Returns the line intensities from two
% intersecting scimat files at a specified resolution. 


%   This function operated on a Slice - by - Slice method: 
%   Scimat1.data and Scimat2.data are 2D datasets 
%   However the data is in 3D space.
%   The user has the choice to specify the Slice and Frame wanted. 
%
%   Function SCIMAT_LINEINTERSECTION() takes two scimat structures
%   along with their coordinates and a specified line resolution. It
%   proceeds to calculate the line of intersection and interpolates the
%   line intensities for both slices. It then returns the line
%   intensities as the outputs. 
%
%   [INTENSITIES1,INTENSITIES2] = SCIMAT_LINEINTERSECTION(SCIMAT1,SCIMAT2,COORDS1,COORDS2,RES)
%       
%   INTENSITIES1,INTENSITIES2 are both the outputs. They are the
%   interpolated line intensities derived from the line between the two 
%   intersecting planes.
%   
%   SCIMAT1,SCIMAT2 (inputs) are both scimat structs
%   
%   COORDS1,COORDS2 (inputs) are the respective cartesian coordinates for
%   scimat1 and scimat2
% 
%   RES is the resolution of the line intensity. It is in mm. A resolution
%   of 1, (RES = 1) accounts for the spacing resolution found in the scimat
%   struct (scimat.axis.spacing). i.e: RES = 1 is equal to a sample point
%   on the line at every 0.0015mm. 
% 
%  
% Authors: Benjamin Villard <b.016434@gmail.com>, Christopher Kelly  <christopher.kelly28@googlemail.com>
% Copyright © 2015 University of Oxford
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



%% Default Checks:

narginchk(2,5)
nargoutchk(0,3)

%Check dimentsion of scimat
chkdims1 = ndims(scimat1.data);
chkdims2 = ndims(scimat2.data);
if chkdims1 > 2 || chkdims2 > 2
    error(' Scimats.data must be 2 Dimensions')
end  

% Check resolution step:

% if isempty (res)
   Res = res; % Set default value of step size top equal 1.40 mm.  
% end



%% First step : Find the direction of the line at the intersection of the planes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculations for scimat1 

% Get corner values to calculate plane vector
S1corner1 = [coords1(1,1,1),coords1(1,1,2),coords1(1,1,3)];
S1corner2 = [coords1(end,1,1),coords1(end,1,2),coords1(end,1,3)];
S1corner3 = [coords1(1,end,1),coords1(1,end,2),coords1(1,end,3)];

S1v1 = S1corner2-S1corner1; % Get plane vector 1
S1v2 = S1corner3-S1corner1; % Get plane vector 2
S1norm = cross(S1v2,S1v1); % Get the orthogonal vector.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculations for scimat2 

% Get corner values to calculate plane vector
S2corner1 = [coords2(1,1,1),coords2(1,1,2),coords2(1,1,3)];
S2corner2 = [coords2(end,1,1),coords2(end,1,2),coords2(end,1,3)];
S2corner3 = [coords2(1,end,1),coords2(1,end,2),coords2(1,end,3)];


S2v1 = S2corner2-S2corner1; % Get plane vector 1
S2v2 = S2corner3-S2corner1; % Get plane vector 2
S2norm = cross(S2v2,S2v1); % Get the orthogonal vector.

%%%%%%%%%%%%%%%%%%%%%%%%Get the direction of the line and make it unitary
% Need to check for amount of slices in volume. 
% Assume there is no in-plane rotation (i.e twisting of plane), Thus Line
% direction will be same throughout slices. LineDir will be computed
% throughout time. 

LineDir = cross(S1norm,S2norm);

if norm(LineDir) < 1e-15
    intensities1 = [];
    intensities2 = [];
else
% make a vector unitary by dividing it with it's norm (length). 
% for F = 1:size(LineDir,4)   
%     LineDir./sqrt(sum(LineDir.^2));
% end

%% Second Step: Calculate a point on the intersecting line

%P = [S1norm(:,:,1,1);S2norm(:,:,1,1);LineDir(:,:,1,1)]/[dot(S1norm(:,:,1,1),S1corner1(:,:,1,1));dot(S2norm(:,:,1,1),S2corner1(:,:,1,1));0]'

% % INPUTS: % % N1 = The first slice's orthogonal vector % N2 = The second
% slice's orthogonal vector % A1 = The first's slice's first corner % A2 =
% The second's slice's first corner % % OUTPUTS: % % P = Point at the
% intersection of the two slices.

% Defaults: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=[0 0 0]; % Initiate point in 3D space
% Point on first plane (take the first corner voxel accross time frames)
PP1 = S1corner1; % Data is 4D (x,y,z accross time)
% Point on second plane (take the first corner voxel accross time frames)
PP2 = S2corner1; % Data is 4D (x,y,z accross time)
% The orthogonal plane vectors are already known:
% S1norm = Orthogonal vector of first planes throughout time (4D) 
% S2norm = Orthogonal vector of second plane throughout time (4D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get the angle between the plane and a point on that line. 
      
d1 = -dot(S1norm,PP1);   %the constants in the Plane 1 equations
d2 = -dot(S2norm, PP2);  %the constants in the Plane 2 equations

% First step is to determine max abs coordinate of cross product
% This allows us to know wether line intersects plane with x = 0,y = 0, or z = 0;
% Next step is to get a point on the intersection line 
% Zero the max coord, and solve for the other two
% There exists three cases which can exists for which the vector containing
% the direction of the line can cross a plane : either the line intersects 
% with x = 0, or y = 0, or z = 0.  The output P gives a point at the
% intersection of the two slices. 

        maxc=find(abs(LineDir)== max(abs(LineDir)));
        
        switch maxc(1)
            case 1                   % intersect with x=0
                P(1)= 0;
                P(2) = (d2*S1norm(3) - d1*S2norm(3))/ LineDir(1);
                P(3) = (d1*S2norm(2) - d2*S1norm(2))/ LineDir(1);
            case 2                    %intersect with y=0
                P(1) = (d1*S2norm(3) - d2*S1norm(3))/ LineDir(2);
                P(2) = 0;
                P(3) = (d2*S1norm(1) - d1*S2norm(1))/ LineDir(2);
            case 3                    %intersect with z=0
                P(1) = (d2*S1norm(2) - d1*S2norm(2))/ LineDir(3);
                P(2) = (d1*S2norm(1) - d2*S1norm(1))/ LineDir(3);
                P(3) = 0;
        end

% P is now a 4D vector containing the x,y,z coordinates for a point on the
% line of intersection between scimat1 and scimat2 across all slices and frames. 

%% Third Step: Calculate coordinates on the intersection lines
% Now we have the direction of the line intersecting both plane, and a
% point on that line. Thus the equation becomes : 
% r =(P(1),P(2),P(3)))+T*(Direction_of_line(1),Direction_of_line(2),Direction_of_line(3))
% Where T in the 'spacing' or 'resolution' of the line (i.e if we want the
% line to be between a discreet interval or to go on infinitely)
% Since the line we get is infinite, we also want to calculate the end points to the line.
% This will be achieved by calculating the points at which the line
% interscect with the sides of the planes (i.e: the vector between their
% corners)


% Defaults: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the longest distance within the image volume is the length of the largest
% diagonal. We compute it in voxel units:
% That is the maximum amount of voxels that can lie on the line. 
lmax = ceil(sqrt(sum(([scimat1.axis.size] - 1).^2))); 
T = 0; % Initiate Resolution; 
% Res = 0.3; % Resolution
lmax = lmax/Res;
lmax = ceil(lmax); % Needed if resolution is greater then 1;
% Preallocate a vector containing the maximum possible amount of points lying
% plane. 
Pmaxp = zeros(lmax,3); % A lmax-by-3 vector to be allocated x, y z coordsflat of line points for NEGATIVE/POSITIVE direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To avoid for loop, we put everything in vector form:

% Vectorize T:
T = repmat((-(lmax-1):Res:lmax-1)',[1,3]);

% Create a matrix of P and LineDir:

P = repmat(P,size(T,1),1);
LineDir = repmat(LineDir,size(T,1),1);

% Calculate line equations:
Pmaxp = P+(T.*(LineDir));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fourth Step: Get intensities from line indices 

% scimat 1
coords1 = reshape(coords1,[numel(coords1(:,:,1)),3])*inv(scimat1.rotmat);

Pmaxp1 = Pmaxp*inv(scimat1.rotmat);

coords1 = reshape(coords1,[size(scimat1.data,1),size(scimat1.data,2),size(scimat1.data,3)*3]);

% scimat 2
coords2 = reshape(coords2,[numel(coords2(:,:,1)),3])*inv(scimat2.rotmat);

Pmaxp2 = Pmaxp*inv(scimat2.rotmat);

coords2 = reshape(coords2,[size(scimat2.data,1),size(scimat2.data,2),size(scimat2.data,3)*3]);

% interpolate scimat1

xv = [ coords1(1,1,1), coords1(1,end,1), coords1(end,end,1),coords1(end,1,1) ];
yv = [ coords1(1,1,2), coords1(1,end,2), coords1(end,end,2),coords1(end,1,2) ];

outOfRange = inpolygon(Pmaxp1(:,1),Pmaxp1(:,2),xv,yv);

Pmaxp1(outOfRange == 0,:) = [];
Pmaxp2(outOfRange == 0,:) = [];

%%%%%%%%%%%%%%%

xv = [ coords2(1,1,1), coords2(1,end,1), coords2(end,end,1),coords2(end,1,1) ];
yv = [ coords2(1,1,2), coords2(1,end,2), coords2(end,end,2),coords2(end,1,2) ];

outOfRange = inpolygon(Pmaxp2(:,1),Pmaxp2(:,2),xv,yv);

Pmaxp1(outOfRange == 0,:) = [];
Pmaxp2(outOfRange == 0,:) = [];

%%%%%%%%%

[coords1(:,:,1),coords1(:,:,2),coords1(:,:,3)] = meshgrid(coords1(1,:,1),coords1(:,1,2),coords1(1,1,3));

intensities1 = interp2( coords1(:,:,1),coords1(:,:,2),...
                        double(scimat1.data),Pmaxp1(:,1),Pmaxp1(:,2));

[coords2(:,:,1),coords2(:,:,2),coords2(:,:,3)] = meshgrid(coords2(1,:,1),coords2(:,1,2),coords2(1,1,3));                    
                    
intensities2 = interp2( coords2(:,:,1),coords2(:,:,2),...
                        double(scimat2.data),Pmaxp2(:,1),Pmaxp2(:,2),'spline'); 
          

 % NEEDED FOR CHRIS' PLOT FUNCTION (validation)
intersectCoords = Pmaxp2*scimat2.rotmat;
end

end
