% SCIMAT  Struct used in Gerardus to store 2D, 3D or 3D+t images and axis
% metainformation.
%
%   Medical image processing usually requires knowledge not only of voxel
%   intensity values, but of voxel size, image orientation and position.
%   Hence, in Gerardus we use a struct, that we call SCIMAT, to enrich
%   images with metainformation. An example of a SCIMAT image:
%
%   scimat = 
%   
%           data:   [415x460x900 single]
%           axis:   [3x1 struct]
%           rotmat: [3x3 double]
%
%   * data:   The image or segmentation.
%   * axis:   Axes metainformation (voxel size, offset and number of voxels).
%   * rotmat: Rotation matrix.
%
% -------------------------------------------------------------------------
%
%   In detail:
%
%   * data:   A 2D to 4D array that contains the image/segmentation voxels.
%
%   * axis:   Axes metainformation. A 1 to 4-vector, depending on the
%             number of dimensions of the image.
%
%     axis(1) ==> rows
%     axis(2) ==> columns
%     axis(3) ==> slices
%     axis(4) ==> time frames
%
%     Each axis element contains the fields:
%
%       size:    Image size.
%       spacing: Voxel size, i.e. spacing between voxel centres. 
%       min:     Coordinates of the "bottom-left" corner of the image, NOT
%                of the first voxel's centre. The first voxel's centre
%                coordinates is offset = min+spacing/2.
%
%   * rotmat: Rotation matrix. A (2, 2) or (3, 3)- matrix for images that
%             are not aligned with the XY- or XYZ-Cartesian axes.
%
%   Note: scimat.axis follows the Matlab convention that in an image
%   (rows, cols) <=> (y, x), but scimat.rotmat is given for (x, y).
%
%   As a historical note, SCIMAT derives from the "scirunnrrd" struct
%   exported by the visualization platform Seg3D, University of Utah
%   Scientific Computing and Imaging (SCI) Institute.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014-2015 University of Oxford
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
