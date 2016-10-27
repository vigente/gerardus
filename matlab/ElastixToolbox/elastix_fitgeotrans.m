function tElx = elastix_fitgeotrans(movingPoints, fixedPoints, transform, sz, spacing, origin)
% ELASTIX_FITGEOTRANS  Transformix struct for image transformation obtained
% from matching two sets of landmarks
%
% The goal of this function is to create a Transformix struct to deform an
% image. The transform maps a set of landmarks on the image to another set
% of landmarks.
%
% TELX = ELASTIX_FITGEOTRANS(MOVINGPOINTS, FIXEDPOINTS, TRANSFORM, SZ, SPACING, ORIGIN)
%
%   TELX is a struct that can be passed to transformix() to deform an
%   image.
%
%   MOVINGPOINTS, FIXEDPOINTS are matrices of the same size, with the
%   coordinates of the target (fixed) and source (moving) landmarks,
%   respectively. Each row contains the coordinates of a point.
%
%   We use Matlab's fitgeotrans() function to find a transformation from
%   the moving to the fixed points. We then invert the transformation
%   matrix, because image transformation goes in the opposite direction as
%   point transformation. Finally, we create the Transformix struct.
%
%   TRANSFORM is a string with the transformation that matches the two sets
%   of landmarks. Currently, we have only implemented 'affine', but it's
%   very easy to add a couple of lines to this function for other
%   transformations supported by Matlab's fitgeotrans ("help fitgeotrans"
%   for options).
%
%   SZ = [cols, rows] is a vector with the size of the output image that
%   transformix will create. Note the Elastix order is different from
%   Matlab's order of size vectors.
%
%   SPACING = [dx, dy] is a vector with the pixel size of the output image.
%
%   ORIGIN = [offx, offy] is a vector with the coordinates of the first
%   pixel in the output image.
%
% See also: transformix, elastix, fitgeotrans.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2016 University of Oxford
% Version: 0.1.0
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
narginchk(6, 6);
nargoutchk(0, 1);

% dimensions of the image and points
D = size(fixedPoints, 2);
if (D ~= size(movingPoints, 2))
    error('movingPoints and fixedPoints must have the same size')
end
if (size(fixedPoints, 1) ~= size(movingPoints, 1))
    error('movingPoints and fixedPoints must have the same size')
end

% compute Matlab transform between two sets of points
tform = fitgeotrans(movingPoints, fixedPoints, transform);

% 'affine'

% basic Elastix transformation struct
tElx.Transform = [];
tElx.NumberOfParameters = NaN;
tElx.TransformParameters = [];
tElx.InitialTransformParametersFileName = 'NoInitialTransform';
tElx.HowToCombineTransforms = 'Compose';
tElx.FixedImageDimension = D;
tElx.MovingImageDimension = D;
tElx.FixedInternalImagePixelType = 'float';
tElx.MovingInternalImagePixelType = 'float';
tElx.Size = sz;
tElx.Index = [0 0];
tElx.Spacing = spacing;
tElx.Origin = origin;
aux = eye(D);
tElx.Direction = aux(:)';
tElx.UseDirectionCosines = 'true';
tElx.CenterOfRotationPoint = [0.0 0.0];
tElx.ResampleInterpolator = 'FinalLinearInterpolator';
tElx.Resampler = 'DefaultResampler';
tElx.DefaultPixelValue = 0;
tElx.ResultImageFormat = 'mha';
tElx.ResultImagePixelType = 'unsigned char';
tElx.CompressResultImage = 'true';

% set transform type
switch (transform)
    case 'affine'
        tElx.Transform = 'AffineTransform';
        tElx.NumberOfParameters = 6;
        
    otherwise
        error(['transform ' transform ' not implemented'])
                   
end

% convert Matlab's transformation matrix to Elastix vector of parameters
% and put it into the basic struct just created
%
% Note that image transformations use the inverse of points
% transformations, to avoid leaving "holes" when interpolating the image
tElx = elastix_affine_matrix2struct(inv(tform.T), tElx);
