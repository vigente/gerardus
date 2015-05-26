function [scimat] = scimat_upsample(scimat,data)

% SCIMAT_UPSAMPLE upsample the scimat image (stored in scimat.data)
% by twice. It first applies a gaussian filter in order to smooth the image
% and the proceeds to sample half of the original image. 
% 
%   DISCLAIMER: Method has NOT been verified for 3D case. 
%    
%   Function SCIMAT_DOWNSAMPLE() takes a scimat struct as an input. It then
%   proceeds to smooth the image in the scimat struct by applying a
%   gaussian filter in either 2 or 3 dimensions, depending on the scimat
%   struct dimensions. Once smoothed, the spacial locations of the voxels
%   is recalcuated at half the sample frequency of the original image. 
%
% [SCIMAT] = SCIMAT_DOWNSAMPLE(SCIMAT)
% 
% SCIMAT (input) and SCIMAT (output) are both Struct used in Gerardus to 
% store 2D, 3D or 3D+t images and axis metainformation. For more
% information see scimat.m
%
% Example: 
%  
% [scimatout] = SCIMAT_DOWNSAMPLE(scimatin)
% 
% scimatin = 
% 
%       axis: [3x1 struct]
%       data: [128x128 double]
%     rotmat: [3x3 double]
%
% scimatout = 
% 
%       axis: [3x1 struct]
%       data: [64x64 double]
%     rotmat: [3x3 double]
%
% Authors: Benjamin Villard <b.016434@gmail.com>,
% Vicente Grau  <vicente.grau@eng.ox.ac.uk>
% Copyright © 2015 University of Oxford
% Version: 0.2.0
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



% Check arguments
narginchk(2, 2);
nargoutchk(0, 1);

% Defaults
orgnlspc(1) = scimat.axis(1).spacing; 
orgnlspc(2) = scimat.axis(2).spacing;
% Get the coordinates from the scimat.
[coords(:,:,1),coords(:,:,2),coords(:,:,3)] = scimat_ndgrid(scimat);
% Obtain the first voxel's spatial position
coordsmin = coords(1,1,:);


% Change size
scimat.axis(1).size = size(data,1);
scimat.axis(2).size = size(data,2);

% Change spacing
scimat.axis(1).spacing = orgnlspc(1)/2;
scimat.axis(2).spacing = orgnlspc(1)/2;

% Need to change the coordinates of first voxel (We use the min
% position of the coordinates and add half of the spacing to get
% the correct scimat struct.
scimat.axis(1).min  = coordsmin(1,1,2) - ([scimat.axis(2).spacing]/2);
scimat.axis(2).min  = coords(1,1,1) - ([scimat.axis(1).spacing]/2);
scimat.axis(3).min  = coords(1,1,3) - ([scimat.axis(3).spacing]/2);

% Update image
scimat.data = data;































