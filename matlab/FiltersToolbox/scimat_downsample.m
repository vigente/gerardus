function [scimat] = scimat_downsample(scimat)
% SCIMAT_DOWNSAMPLE downsamples the scimat image (stored in scimat.data)
% by half. It first applies a gaussian filter in order to smooth the image
% and the proceeds to sample half of the original image. 
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
narginchk(1, 1);
nargoutchk(0, 1);


if(ndims(scimat) == 3)
    
    % Create Gaussian filter
    gfilt = zeros(1,1,5);   %Three-dim kernel
    gfilt(1,1,:) = [1 4 6 4 1]./16;
    
    % Apply filter to image in all three directions
    scimat = convn(convn(convn(scimat,shiftdim(gfilt,2),'same'),shiftdim(gfilt,1),'same'),gfilt,'same');
    
    % Recalculate image
    scimatout = scimat(1:2:size(scimat,1),1:2:size(scimat,2),1:2:size(scimat,3));
     
    % Update scimat
    scimat.data = scimatout;
    
    % Recalculate size of image
    scimat.axis(1).size = size(scimat.data,1);
    scimat.axis(2).size = size(scimat.data,2);
    scimat.axis(3).size = size(scimat.data,3);

elseif(ndims(scimat.data) == 2)
    
    % Create Gaussian filter
    gfilt=[1 4 6 4 1]./16;
    
    % Apply filter to image in both directions
    scimat.data = convn(scimat.data,gfilt,'same');
    scimat.data = convn(scimat.data,gfilt','same');        
    
    % Recalculate size of image
    size_out = [floor((scimat.axis(1).size)/2) floor((scimat.axis(2).size)/2)];
    
    % Update scimat
    scimat.axis(1).size = size_out(1);
    scimat.axis(2).size = size_out(2);
    
    % Recalculate image
    [x,y] = ndgrid(1.5:2:(size_out(1)*2)-.5, 1.5:2:(size_out(2)*2)-.5);
    scimatout=interpn(scimat.data,x,y);
    
    % Update scimat
    scimat.data = scimatout;
    
end
end