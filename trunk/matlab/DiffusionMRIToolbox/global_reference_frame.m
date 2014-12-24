function [ rad, circ, long ] = global_reference_frame(LV_cavity, orientation, dims)
%GLOBAL_REFERENCE_FRAME Returns unit vectors for the global reference frame
%   Fits a straight line through the center of the left ventricle cavity,
%   and then computes the radial, circumferential, and longitudinal vectors
%   based on this line.
%
%
%   Inputs:
%       LV_CAVITY is a 3D binary segmentation of the left ventricle cavity.
%           The long axis should be more or less aligned with the 3rd dimension
%           of the image (z axis).
%       ORIENTATION can be: 'ApexFirstAxial' 
%                           'ApexLastAxial' (default)
%                           'ApexFirstSagittal'
%                           'ApexLastSagittal'
%                           'ApexFirstCoronal'
%                           'ApexLastCoronal'
%       DIMS is the voxel size of the data (default [1 1 1])
%
%   Outputs:
%       RAD is the radial unit vectors (along the 4th dimension)
%       CIRC is the circumferential unit vectors
%       LONG is the longitudinal unit vectors


% Author: Darryl McClymont <darryl.mcclymont@gmail.com>
% Copyright � 2014 University of Oxford
% Version: 0.1.3
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
narginchk(1,3);
nargoutchk(0, 3);

if nargin < 2
    orientation = 'ApexLastAxial';
end
if nargin < 3
    dims = [1 1 1];
end

% check orientation
needs_permuting = 0;
if strcmp(orientation, 'ApexFirstSagittal') || ...
    strcmp(orientation, 'ApexLastSagittal')
    
    needs_permuting = 1;
    
    perm_mat = [2 3 1];
    inv_perm_mat = [3 1 2];
    
    dims = dims(perm_mat);
    LV_cavity = permute(LV_cavity, perm_mat);

elseif strcmp(orientation, 'ApexFirstCoronal') || ...
    strcmp(orientation, 'ApexLastCoronal')
    
    needs_permuting = 1;
    
    perm_mat = [1 3 2];
    inv_perm_mat = [1 3 2];
    
    dims = dims(perm_mat);
    LV_cavity = permute(LV_cavity, perm_mat);

end

% when apex is last, the longitudinal vector is flipped
long_multiplier = 1;
if strcmp(orientation, 'ApexLastSagittal') || ...
    strcmp(orientation, 'ApexLastCoronal') || ...
    strcmp(orientation, 'ApexLastAxial') 

    long_multiplier = -1;
end


sz = size(LV_cavity);

row_vector = zeros(1,sz(3));
col_vector = zeros(1,sz(3));
slice_vector = 1:sz(3);

% run though all slices, and fit a line through the center of the LV cavity
for i = 1:size(LV_cavity,3)
   C = LV_cavity(:,:,i);
   if sum(C(:)) == 0
       continue
   end
   
   cr = getfield(regionprops(C, 'Centroid'), 'Centroid');
   
   row_vector(i) = cr(2);
   col_vector(i) = cr(1);
      
end

% Only keep those with non-zero voxel counts
fitrow = row_vector(row_vector ~= 0) * dims(1);
fitcol = col_vector(row_vector ~= 0) * dims(2);
fitslice = slice_vector(row_vector ~= 0) * dims(3);


% These lines copied from Valentina's ScriptAssignGlobalCylindricalSysToPhantom
% fit a straight line through the centroids
n=1;%linear
slope_xz = polyfit(fitslice,fitrow,n);
slope_yz =polyfit(fitslice,fitcol,n);

% origin
Origin = [slope_xz(2), slope_yz(2), 0];
O = repmat(permute(Origin, [1 3 4 2]), [sz, 1]);

% longitude vector
l = [slope_xz(1), slope_yz(1), 1];
l = long_multiplier * l / norm(l);
long = repmat(permute(l, [1 3 4 2]), [sz, 1]);

% voxel coordinates
[X, Y, Z] = ndgrid(1:sz(1), 1:sz(2), 1:sz(3));
V = cat(4, X, Y, Z);

% use the dot product to find the position along the longitudinal vector
long_mag = dot(V - O, long, 4);
vox_nearest_point = O + long .* repmat(long_mag, [1 1 1 3]);

%         ^ longitudinal vector
%        /
%       /
%      . voxel nearest point
%     / \ 
%    /   \ radial vector
%   /     . voxel coordinates
%  /
% /
%. origin

% radial vector
rad = V - vox_nearest_point;
rad = rad ./ repmat(sqrt(sum(rad.^2,4)), [1 1 1 3]);

% circumferential vector = long X rad
circ = cross(long, rad, 4);
circ = circ ./ repmat(sqrt(sum(circ.^2,4)), [1 1 1 3]);


if needs_permuting
    rad = permute(rad(:,:,:,inv_perm_mat), [inv_perm_mat, 4]);
    circ = permute(circ(:,:,:,inv_perm_mat), [inv_perm_mat, 4]);
    long = permute(long(:,:,:,inv_perm_mat), [inv_perm_mat, 4]);
end

end

