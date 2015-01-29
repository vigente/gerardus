function scimat = dcm2scimat(dcm,dcminfo)
% DCM2SCIMAT  converts dicom datastructure to scimat format.
%
% SCIMAT = DCM2SCIMAT(DCM, DCMINFO)
%
%   DCM is a 2D/3D, 4-D array with the DICOM images. The format is 
%   DCM(row, column, slice, time frame).
%
%   DCMINFO is a matrix of structs. DCMINFO(slice, time frame) contains the
%   metainformation for the corresponding slice and time frame in the data
%   set.
%
%   SCIMAT is a Struct used in Gerardus to store 3D images and their
%   metadata. Please see scimat.m for more info. 
%   
%   SCIMAT is a struct wich returns 3 fields:    
%   
%           data:   An array containing the dicom images, in a volume format. Data
%                   can be 2D, 3D, or 4D (3D + time). 
%           axis:   A 3-vector with the following fields:
%                   
%                   size: Image size.
%                    
%                   spacing: Voxel size, i.e. spacing between voxel centres. 
%                    
%                   min: Coordinates of the "bottom-left" corner of the image, NOT
%                       of the first voxel's centre. The first voxel's centre
% 
%           rotmat: A (3, 3)- matrix that gives the image orientation as a rotation
%                   of the x, y, z Cartesian axes.
%   
%   The SCIMAT struct follows the Matlab convention that image arrays are
%   sorted as (y, x, z)
%
%     axis(1) ==> rows ==> y axis
%     axis(2) ==> columns ==> x axis
%     axis(3) ==> slices
% 
% See also: scimat.m.
%
% Authors:Benjamin Villard <b.016434@gmail.com>
% Copyright © 2014-2015 University of Oxford
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
% http://code.google.com/p/gerardus/
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
%   


% check arguments
narginchk(1, 2);
nargoutchk(0, 1);

% defaults
if (nargin < 1 || isempty(dcm))
    dcm = [];
end
if (nargin < 2 || isempty(dcminfo))
    dcminfo = [];
end


% data volume (Dicom Image)
scimat.data = dcm;

scimat.axis(1).size = size(dcm,1);% Resolution in mm
scimat.axis(2).size = size(dcm,2);
scimat.axis(3).size = size(dcm,3); 

% Check if Slices are 3D. NOTE: the information concerning the spacing of
% the pixel/voxel was taken from the first file info. This is because we
% assume a consistency in the imaging acquisition. The resolution of the
% image should not change in space and time. 
if isempty(dcminfo(1,1).SliceThickness)
    vol = 1; % Slice Thickness is NOT empty --> Voxel (3D)
else
    vol = 2; % Slice Thickness is empty --> Pixel (2D)
end

switch vol
    case 1
        scimat.axis(1).spacing = dcminfo(1,1).PixelSpacing(1);% Resolution in mm
        scimat.axis(2).spacing = dcminfo(1,1).PixelSpacing(2);
        scimat.axis(3).spacing = 1; 
        
    case 2
        scimat.axis(1).spacing = dcminfo(1,1).PixelSpacing(1);% Resolution in mm
        scimat.axis(2).spacing = dcminfo(1,1).PixelSpacing(2);
        scimat.axis(3).spacing = dcminfo(1,1).SliceThickness; 
        
end

% left edge of first voxel. This metric is different for each Slice, but consistent
% throughout time. Hence min will be a 1xS array of the x,y,z positions of
% the slices. Also we assume notation is the following: data(rows,columns,slice,time). 
% As slice position is equal throughout time, T = 1;  

for S = 1:size(scimat.data,3)
    scimat.axis(1).min(S) = dcminfo(S,1).ImagePositionPatient(2); % min Y position of the left edge of the first voxel 
    scimat.axis(2).min(S) = dcminfo(S,1).ImagePositionPatient(1); % min X position of the left edge of the first voxel 
    scimat.axis(3).min(S) = dcminfo(S,1).ImagePositionPatient(3); % min Z position of the left edge of the first voxel 
end
% Rotation Matrix. Same throughout space and time. 
    scimat.rotmat = [ dcminfo(1,1).ImageOrientationPatient(1:3),...
                                                    dcminfo(1,1).ImageOrientationPatient(4:6),...
                                                    cross(dcminfo(1,1).ImageOrientationPatient(1:3), dcminfo(1,1).ImageOrientationPatient(4:6)) ];

scimat.axis = scimat.axis'; 
end 

