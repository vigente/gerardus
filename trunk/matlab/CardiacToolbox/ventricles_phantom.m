function [im,lv,rv,myo] = ventricles_phantom(m,v)
% VENTRICLES_PHANTOM Create phantom image of simplified heart ventricles
%
% [IM,LV,RV,MYO] = ventricles(M,V)
%
%   IM is a 3D image that consists of a heart phantom depicting the left
%   ventricle, right ventricle and myocardium
%
%   LV is a 3D image that consists of the LV cavity phantom
%
%   RV is a 3D image that consists of the RV cavity phantom
%
%   MYO is a 3D image that consists of a myocardium phantom
%
%   M is the mean Gaussian white noise added to the final heart phantom
%   generated (the default value is zero)
%
%   V is the the variance of the Gaussian white noise added to the final
%   heart phantom generated (the default value is zero)

% Author: Christopher Kelly <christopher.kelly28@gmail.com>
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

narginchk(0,2);
nargoutchk(0,4);

% default values

if (nargin == 1)
    error('You must either specify all 2 input parameters or none at all')
elseif (nargin < 2)
    m      = 0;
    v      = 0;
end

% TODO: we hard-code these values. In the future, we may want to have these
% as input arguments, but currently changing these values produces odd
% looking results
lvendo = 300;
lvepi  = 400;
rvendo = 300;
rvepi  = 400;

%% create image space

xRange = -(128/2) + 0.5 : (128/2) - 0.5;
yRange = xRange;
zRange = xRange;

[ x, y, z ] = ndgrid( xRange, yRange, zRange );

im = zeros(size(x));

lv = im;
rv = im;
myo = im;
heart = im;

%% define ellipsoids shapes

radius1 = ( (x).^2 + (y).^2 + (z/2).^2 );
radius2 = ( (x).^2 + (y/1.7).^2 + (z/2).^2 );

%% define voxels in left ventricle wall

% define LV ellpsoid edge thickness
ellipsoid1Outer = find(radius1 < lvepi );
ellipsoid1Inner = find(radius1 < lvendo );

% define LV ellipsoid
ellipsoid1 = ismember(ellipsoid1Outer,ellipsoid1Inner);
ellipsoid1 = ellipsoid1Outer(ellipsoid1 == 0);

[ellipsoid1X, ellipsoid1Y,ellipsoid1Z] = ind2sub(size(im),ellipsoid1);

ellipsoid1X(ellipsoid1Z > 70) = [];
ellipsoid1Y(ellipsoid1Z > 70) = [];
ellipsoid1Z(ellipsoid1Z > 70) = [];

% find voxel numbers that will make up LV wall
LV_wall = sub2ind(size(im),ellipsoid1X,ellipsoid1Y,ellipsoid1Z);

%% define voxels in right ventricle wall

% define RV ellpsoid edge thickness
ellipsoid2Outer = find(radius2 < rvepi );
ellipsoid2Inner = find(radius2 < rvendo );

% define RV ellipsoid
ellipsoid2 = ismember(ellipsoid2Outer, ellipsoid2Inner);
ellipsoid2 = ellipsoid2Outer(ellipsoid2 == 0);                                	% find voxel numbers that will make up RV wall

[ellipsoid2X, ellipsoid2Y,ellipsoid2Z] = ind2sub(size(im),ellipsoid2);

ellipsoid2X(ellipsoid2Z > 70) = [];
ellipsoid2Y(ellipsoid2Z > 70) = [];
ellipsoid2Z(ellipsoid2Z > 70) = [];

ellipsoid2X(ellipsoid2Y > (size(im,2)/2)) = [];
ellipsoid2Z(ellipsoid2Y > (size(im,2)/2)) = [];
ellipsoid2Y(ellipsoid2Y > (size(im,2)/2)) = [];

RV_wall = sub2ind(size(im),ellipsoid2X,ellipsoid2Y,ellipsoid2Z);

%% define left ventricle

lv(radius1 < lvendo) = 1;

LVellipsoid = find(lv == 1);
[LVellipsoidX, LVellipsoidY, LVellipsoidZ] = ind2sub(size(im),LVellipsoid);

LVellipsoidX(LVellipsoidZ < 70) = [];
LVellipsoidY(LVellipsoidZ < 70) = [];
LVellipsoidZ(LVellipsoidZ < 70) = [];

LVellipsoid = sub2ind(size(im),LVellipsoidX, LVellipsoidY, LVellipsoidZ);
lv(LVellipsoid) = 0;

%% define right ventricle

rv(radius2 < rvendo) = 1;

RVellipsoid = find(rv == 1);
[RVellipsoidX, RVellipsoidY, RVellipsoidZ] = ind2sub(size(im),RVellipsoid);

RVellipsoidX( (RVellipsoidZ < 70) ) = [];
RVellipsoidY( (RVellipsoidZ < 70) ) = [];
RVellipsoidZ( (RVellipsoidZ < 70) ) = [];

RVellipsoid = sub2ind(size(im),RVellipsoidX, RVellipsoidY, RVellipsoidZ);
rv(RVellipsoid) = 0;

RVellipsoid = find(rv == 1);
[RVellipsoidX, RVellipsoidY, RVellipsoidZ] = ind2sub(size(im),RVellipsoid);

RVellipsoidX(RVellipsoidY < (size(im,2)/2)) = [];
RVellipsoidZ(RVellipsoidY < (size(im,2)/2)) = [];
RVellipsoidY(RVellipsoidY < (size(im,2)/2)) = [];

RVellipsoid = sub2ind(size(im),RVellipsoidX, RVellipsoidY, RVellipsoidZ);
rv(RVellipsoid) = 0;

rv([LV_wall;find(lv == 1)]) = 0;

%% define myocardium

myo([LV_wall;RV_wall]) = 1;

%% define heart

heart( lv == 1 ) = 1;
heart( rv == 1 ) = 2;
heart( myo == 1 ) = 3;

%% add noise

heart = (heart - min(heart(:)))/(max(max(max(heart)))-min(min(min(heart))));

for i = 1:size(im,3)
    heart(:,:,i) = imnoise(heart(:,:,i),'Gaussian',m,v);
    lv(:,:,i) = imnoise(lv(:,:,i),'Gaussian',m,v);
    rv(:,:,i) = imnoise(rv(:,:,i),'Gaussian',m,v);
    myo(:,:,i) = imnoise(myo(:,:,i),'Gaussian',m,v);
end

im = heart;





