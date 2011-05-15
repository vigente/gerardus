function stats = scinrrd_regionprops(nrrdl, varargin)
% SCINRRD_REGIONPROPS  Measure properties of image regions on each slice
%
%   This function basically applies Matlab's regionprops() function to each
%   slice of an SCI NRRD volume.
%
% STATS = SCINRRD_REGIONPROPS(NRRDL, PROPERTIES)
%
%   NRRDL is the struct with the segmentation. As in Matlab's function
%   regionprops(), each positive integer value of NRRDL corresponds to a
%   different region.
%
% STATS = SCINRRD_REGIONPROPS(NRRDL, NRRDI, PROPERTIES)
%
%   NRRDI is a struct with a greyscale volume instead of a segmentation.
%   Both NRRDL and NRRDI must have the same size.
%
%   From the regionprops()'s help: 
%
%   PROPERTIES can be a comma-separated list of strings, a cell array
%   containing strings, the string 'all', or the string 'basic'. The set of
%   valid measurement strings includes:
% 
%   Shape measurements
% 
%     'Area'              'EulerNumber'       'Orientation'               
%     'BoundingBox'       'Extent'            'Perimeter'          
%     'Centroid'          'Extrema'           'PixelIdxList' 
%     'ConvexArea'        'FilledArea'        'PixelList'
%     'ConvexHull'        'FilledImage'       'Solidity' 
%     'ConvexImage'       'Image'             'SubarrayIdx'            
%     'Eccentricity'      'MajorAxisLength' 
%     'EquivDiameter'     'MinorAxisLength'                   
%      
%   Pixel value measurements (requires grayscale image as an input)
% 
%     'MaxIntensity'
%     'MeanIntensity'
%     'MinIntensity'
%     'PixelValues'
%     'WeightedCentroid'
% 
%   Property strings are case insensitive and can be abbreviated.
%
%   If PROPERTIES is the string 'all', REGIONPROPS returns all of the
%   Shape measurements. If called with a grayscale image, REGIONPROPS also
%   returns Pixel value measurements. If PROPERTIES is not specified or if
%   it is the string 'basic', these measurements are computed: 'Area',
%   'Centroid', and 'BoundingBox'.
%
%
%   Note on SCI NRRD: Software applications developed at the University of
%   Utah Scientific Computing and Imaging (SCI) Institute, e.g. Seg3D,
%   internally use NRRD volumes to store medical data.
%
%   When label volumes (segmentation masks) are saved to a Matlab file
%   (.mat), they use a struct called "scirunnrrd" to store all the NRRD
%   information:
%
%   >>  scirunnrrd
%
%   scirunnrrd = 
%
%          data: [4-D uint8]
%          axis: [4x1 struct]
%      property: []

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010 University of Oxford
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
error( nargchk( 1, Inf, nargin, 'struct' ) );
error( nargoutchk( 0, 1, nargout, 'struct' ) );

% parse input
if (nargin > 1)
    if (isstruct(varargin{1})) 
        if (nargin == 2) % SCINRRD_REGIONPROPS(NRRDL, NRRDI)
            nrrdi = varargin{1};
            args = [];
        else % SCINRRD_REGIONPROPS(NRRDL, NRRDI, PROPERTIES)
            nrrdi = varargin{1};
            args = varargin(2:end);
        end
    else % SCINRRD_REGIONPROPS(NRRDL, PROPERTIES)
        nrrdi = [];
        args = varargin(:);
    end
end

% if we have two volumes, they must be the same size
if ~isempty(nrrdi)
    if ~isequal(size(nrrdl.data), size(nrrdi.data))
        error('NRRDL and NRRDI data must have the same size')
    end
end

% init output
stats = cell(size(nrrdl.data,3), 1);

% go slice by slice
for I = 1:size(nrrdl.data,3)
    if ~isempty(nrrdi)
        stats{I} = regionprops(nrrdl.data(:,:,I), nrrdi.data(:,:,I), args);
    else
        stats{I} = regionprops(nrrdl.data(:,:,I), args);
    end
end
