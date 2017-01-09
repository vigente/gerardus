function [ h ] = quiver_image( IM, VEC, orient, ix , VEC2)
%QUIVER_IMAGE Displays a 2D slice of a 3D image with a quiver plot
%   overlaid on top.
%
%   Input: 
%       IM: 3D image (row col slice)
%       VEC: 4D vector field (row col slice vector), where the vector
%           dimension is of length 3 (x y z)
%       ORIENT: which plane do you want to display? A string of 'r', 'c', 
%           or 's'
%       IX: the index of the row/column/slice to display (integer)
%       VEC2: if you want to have a second vector field on top, put it here
%           (but it gets a bit crowded)
%
%   Output: 
%       H: figure handle

% Author: Darryl McClymont  <darryl.mcclymont@gmail.com>
% Copyright © 2015 University of Oxford
% Version: 0.1.0
%
%University of Oxford means the Chancellor, Masters and Scholars of
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


if nargin < 3
    orient = 'slice';
end

switch orient
    case{'row','r'}

        dim_to_collapse = 1;

    case{'col', 'c', 'column'}

        dim_to_collapse = 2;

    case{'s', 'slice'}

       dim_to_collapse = 3;

    otherwise         
        error('Orient should be a string, equal to row, r, col, c, column, s, or slice')
end


if nargin < 4

    ix = ceil(size(IM,dim_to_collapse)/2);
    
end
  

VEC(isinf(VEC)) = 0;

null_mask = double(max(abs(VEC), [], 4) == 0);
null_mask(null_mask == 1) = NaN;
null_mask(null_mask == 0) = 1;
VEC = bsxfun(@times, VEC, null_mask);

% reshape the image and vector field such that the collapseable dimension
% is at the end

dim_order1 = [setdiff(1:3, dim_to_collapse), dim_to_collapse];
IM = permute(IM, dim_order1);
IM = IM(:,:,ix);

dim_order = [setdiff(1:4, dim_to_collapse), dim_to_collapse];
VEC = permute(VEC, dim_order);
VEC = VEC(:,:,:,ix);

if nargin > 4
    VEC2 = permute(VEC2, dim_order);
    VEC2 = VEC2(:,:,:,ix);
end

xstep = size(IM,1)/size(VEC,1);
ystep = size(IM,2)/size(VEC,2);

[X, Y] = ndgrid(xstep/2:xstep:size(IM,1), ystep/2:ystep:size(IM,2));

AxesHandle = gca;
if ~isempty(get(AxesHandle, 'Children')) % if the current axis is not empty
    figure;
end
h = imagesc(IM);
set(h,'alphadata',~isnan(VEC(:,:,dim_order1(1))));
set(gca, 'Color', [0 0 0])



colormap jet
hold on; 

% set the top left voxel to something high, so that the others are scaled
% better
VEC(1,1,dim_order1(1)) = max(VEC(:))*sqrt(2);
VEC(1,1,dim_order1(2)) = max(VEC(:))*sqrt(2);

quiver(Y+ystep/2,X+xstep/2, VEC(:,:,dim_order1(2)), VEC(:,:,dim_order1(1)), 'k', 'LineWidth', 1.2);
quiver(Y+ystep/2,X+xstep/2, -VEC(:,:,dim_order1(2)), -VEC(:,:,dim_order1(1)), 'k', 'LineWidth', 1.2);

if nargin > 4
    quiver(Y+ystep/2,X+xstep/2, VEC2(:,:,dim_order1(2)), VEC2(:,:,dim_order1(1)), 'w', 'LineWidth', 1.2);
    quiver(Y+ystep/2,X+xstep/2, -VEC2(:,:,dim_order1(2)), -VEC2(:,:,dim_order1(1)), 'w', 'LineWidth', 1.2);
end

axis image


end

