function scimat = dcm2scimat(dcm, dcminfo,'mode')
% DCM2SCIMAT  Load DICOM file or convert DICOM image data into SCIMAT
% format.
%
% SCIMAT = DCM2SCIMAT(FILENAME)
%
%   FILENAME is a string with the path and name of a DICOM file.
%
%   SCIMAT is a struct used in Gerardus to store images and their metadata
%   (see scimat.m for more details).
%
%   Note 1: DICOM uses mm and ms as spatial and temporal units, but we
%   convert them to metres and seconds, respectively.
%
%   Note 2: DICOM direction cosines define a left-multiplying rotation
%   matrix. We transpose it to make it right-multiplying, because in Matlab
%   points are usually given as row vectors.
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
%   MODE is a string. The string can take on two options:
%                         'volume': Creates a uniform volume 
%                         'slice': Creates multiple scimat files of the
%                         form: 
%                                 scimat(s).data(:,:,1,:)
%                                 scimat(s).axis(:,:,1,:)
%                                 scimat(s).rotmat
%                                     
%    Note: If the volume is non uniform, the program automatically
%    executes as a slice method. 
% 
% See also: scimat.m.

% Authors: Benjamin Villard <b.016434@gmail.com>,
% Ramon Casero <rcasero@gmail.com>
% Copyright © 2014-2015 University of Oxford
% Version: 0.4.1
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

% check arguments
narginchk(1, 3);
nargoutchk(0, 1);

% if the user provides only one input argument, we assume that it is a
% filename to load
if (nargin == 1)
    filename = dcm;
    dcm = dicomread(filename);
    dcminfo = dicominfo(filename);
end

% check inputs
if (size(dcm, 3) ~= size(dcminfo, 1))
    error('Number of slices in DCM must be the same as number of rows in DCMINFO')
end
if (size(dcm, 4) ~= size(dcminfo, 2))
    error('Number of frames in DCM must be the same as number of columns in DCMINFO')
end


% number of rows, columns, slices and frames
NR = size(dcm, 1);
NC = size(dcm, 2);
NS = size(dcm, 3);
NF = size(dcm, 4);

%size of original position vector (should always be 3 for x,y,z)
SV = size(dcminfo(1, 1).ImagePositionPatient, 1); 

% Offset
offset = dcminfo(1, 1).ImagePositionPatient;

if strcompi(varargin{3},'volume')
   
%% Check for uniform spacing
% First make sure slices are ordered: 
[~,locationOrder] = sort([dcminfo(:,1).SliceLocation]);
dcminfo = dcminfo(locationOrder,:);
dcm = dcm(:,:,locationOrder,:);

% Check for uniformity in spacing 

pos = zeros(SV, NS); % Pre-allocate space 

% Store x,y,z position of upper left corner of image voxel (center of voxel) 
% for all slices

for I = 1:NS
    pos(:, I) = dcminfo(I, 1).ImagePositionPatient(:); 
end

d = diff(pos, 1, 2); % get the difference between each point (x,y,z) in 
                     % space from one slice to the next. 

chkd = sqrt(sum(d.^2, 1)); % Get the norm between all the slices. (Norm of vectors)

% Check if some slices are retaken, and if so, remove them
idx = find(chkd ==0);
chkd(idx) = [];
dcm(:,:,idx,:) = [];
dcminfo(:,:,idx,:) = [];

% Check if same spacing between slices. % If No error is given then code
% assumes volume is uniform and continues. 
for R = 1:size(chkd,2)
   
    if round(chkd(R)) ~= round(mode(chkd))
        error('Space between slices is non-uniform!')
    end
    
end

end   

%% Build scimat. 
% At this point volume should be uniform. 

%% rotation matrix [Builds Rotmat, and also checks if rotation matrix is orthonormal]. 
rotmat = [ ...
    dcminfo(1,1).ImageOrientationPatient(1:3),... % 1st column before transpose
    dcminfo(1,1).ImageOrientationPatient(4:6),... % 2nd column
    cross(dcminfo(1,1).ImageOrientationPatient(1:3), ... % 3rd column
    dcminfo(1,1).ImageOrientationPatient(4:6)) ...
    ]';

% if (NS == 1) % one slice
     
% else % multiple slices    
%     % 3rd column of rotation matrix (before transposing) scaled by the voxel
%     % size
%     rot3 = (dcminfo(end, 1).ImagePositionPatient ...
%          - dcminfo(1, 1).ImagePositionPatient)'/(NS-1) * 1e-3;
%     
%     % voxel size in the slice dimension
%     rot3 ./ rotmat(3, :)  
% end

%% image metadata
scimat.axis(1).size = size(dcm,1);
scimat.axis(2).size = size(dcm,2);
scimat.axis(3).size = size(dcm,3);
if (ndims(dcm) > 3)
    scimat.axis(4).size = size(dcm,4); 
else
    scimat.axis(4).size = 1;
end

% Check if Slices are 3D. NOTE: the information concerning the spacing of
% the pixel/voxel was taken from the first file info. This is because we
% assume a consistency in the imaging acquisition. The resolution of the
% image should not change in space and time. 
if isempty(dcminfo(1,1).SliceLocation) % Slice Location is empty --> Pixel (2D)
    scimat.axis(1).spacing = dcminfo(1,1).PixelSpacing(1) * 1e-3; % units: m
    scimat.axis(2).spacing = dcminfo(1,1).PixelSpacing(2) * 1e-3; % units: m
    scimat.axis(3).spacing = 1;
    % check if data is temporal
    if (ndims(scimat.data) > 3)
        scimat.axis(4).spacing = dcminfo(1,2).TriggerTime * 1e-3; % units: seconds
    else
        scimat.axis(4).spacing = 1;
    end
else
    % Slice Thickness is NOT empty --> Voxel (3D)
    scimat.axis(1).spacing = dcminfo(1,1).PixelSpacing(1) * 1e-3; % units: m
    scimat.axis(2).spacing = dcminfo(1,1).PixelSpacing(2) * 1e-3; % units: m
    %Check if multiple slices or whether data is just one slice but
    %temporal.
    for S = 1:size(dcminfo,1)
        if size(dcminfo,1) == 1 % So that NaN errors don't come up due to the mean(diff())
            scimat.axis(3).spacing = abs(mean([dcminfo(:,1).SliceLocation])) * 1e-3; % units: m 
        else
            scimat.axis(3).spacing = abs(mean(diff([dcminfo(:,1).SliceLocation]))) * 1e-3; % units: m
        end
    end
    %Check if data is temporal
    if (ndims(dcm) > 3)
        scimat.axis(4).spacing = dcminfo(1,2).TriggerTime * 1e-3; % units: seconds
    else
        scimat.axis(4).spacing = 1;
    end
end

% left edge of first voxel. This metric is different for each Slice, but consistent
% throughout time. Hence min will be a 1xS array of the x,y,z positions of
% the slices. Also we assume notation is the following: data(rows,columns,slice,time). 
% As slice position is equal throughout time, T = 1;  

scimat.axis(1).min = offset(2) * 1e-3 - scimat.axis(1).spacing/2; % min Y position of the left edge of the first voxel
scimat.axis(2).min = offset(1) * 1e-3 - scimat.axis(2).spacing/2; % min X position of the left edge of the first voxel
scimat.axis(3).min = offset(3) * 1e-3 - scimat.axis(3).spacing/2; % min Z position of the left edge of the first voxel
if (ndims(dcm) > 3)
    scimat.axis(4).min = dcminfo(1,1).TriggerTime  * 1e-3;% min of Trigger time is 0.
else
    scimat.axis(4).min = 0;% min of Trigger time is 0.
end




%% create the scimat struct

% image volume (Dicom Image)
scimat.data = dcm;

% Rotation Matrix. Same throughout space and time. 
scimat.rotmat = rotmat;

% Axis 
scimat.axis = scimat.axis'; 

