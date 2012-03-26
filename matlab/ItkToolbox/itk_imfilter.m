function im = itk_imfilter(~, ~)
% ITK_IMFILTER: Run ITK filter on a 2D or 3D image
%
% This MEX function is a multiple-purpose wrapper to be able to run
% all ITK filters that inherit from itk::ImageToImageFilter on a
% Matlab 2D image or 3D image volume.
%
% B = ITK_IMFILTER(TYPE, A, [FILTER PARAMETERS])
%
%   TYPE is a string with the filter we want to run. See below for a whole
%   list of options.
%
%   A is a 2D matrix or 3D volume with the image or
%   segmentation. Currently, A can be of any of the following
%   Matlab classes:
%
%     boolean
%     double
%     single
%     int8
%     uint8
%     int16
%     uint16
%     int32
%     int64
%
%   A can also be a SCI NRRD struct, A = nrrd, with the following fields:
%
%     nrrd.data: 2D or 3D array with the image or segmentation, as above
%     nrrd.axis: 3x1 struct array with fields:
%       nnrd.axis.size:    number of voxels in the image
%       nnrd.axis.spacing: voxel size, image resolution
%       nnrd.axis.min:     real world coordinates of image origin
%       nnrd.axis.max:     ignored
%       nnrd.axis.center:  ignored
%       nnrd.axis.label:   ignored
%       nnrd.axis.unit:    ignored
%
%   (An SCI NRRD struct is the output of Matlab's function scinrrd_load(),
%   also available from Gerardus.)
%
%   [FILTER PARAMETERS] is an optional list of parameters, specific for
%   each filter. See below for details.
%
%   B has the same size as the image in A, and contains the filtered image
%   or segmentation mask. It's type depends on the type of A and the filter
%   used, and it's computed automatically.
%
%
% Supported filters:
% -------------------------------------------------------------------------
%
% B = ITK_IMFILTER('skel', A)
%
%   (itk::BinaryThinningImageFilter3D) Skeletonize a binary mask
%
%   B has the same size and class as A
%
% [B, NV] = ITK_IMFILTER('dandist', A)
%
%   (itk::DanielssonDistanceMapImageFilter) Compute unsigned distance map
%   for a binary mask. Distance values are given in voxel coordinates
%
%   B has the same size as A. B has a type large enough to store the
%   maximum distance in the image. The largest available type is double. If
%   this is not enough, a warning message is displayed, and double is used
%   as the output type.
%
%   NV has the same size as A. Each element has the index of the closest
%   foreground voxel. For example, NV(4) = 7 means that voxel 4 is the
%   closest foreground voxel to voxel 7.
%
% B = ITK_IMFILTER('maudist', A)
%
%   (itk::SignedMaurerDistanceMapImageFilter) Compute signed distance map
%   for a binary mask. Distance values are given in real world coordinates,
%   if the input image is given as an NRRD struct, or in voxel units, if
%   the input image is a normal array. 
%
%   The output type is always double.
%
% B = ITK_IMFILTER('bwdilate', A, RADIUS, FOREGROUND)
% B = ITK_IMFILTER('bwerode', A, RADIUS, FOREGROUND)
%
%   (itk::BinaryDilateImageFilter). Binary dilation. The structuring
%   element is a ball.
%   (itk::BinaryErodeImageFilter). Binary erosion. The structuring
%   element is a ball.
%
%   RADIUS is a scalar with the radius of the ball in voxel units. If a
%   non-integer number is provided, then floor(RADIUS) is used. By default,
%   RADIUS = 0 and no dilation is performed.
%
%   FOREGROUND is a scalar. Voxels with that value will be the only ones
%   dilated. By default, FOREGROUND is the maximum value allowed for the
%   type, e.g. FOREGROUND=255 if the image is uint8. This is the default in
%   ITK, so we respect it.
%
% B = ITK_IMFILTER('advess', A, SIGMAMIN, SIGMAMAX, NUMSIGMASTEPS, NUMITERATIONS,
%                  WSTRENGTH, SENSITIVITY, TIMESTEP, EPSILON)
%
%   (itk::AnisotropicDiffusionVesselEnhancementImageFilter)
%   Anisotropic difussion vessel enhancement.
%
%   Enquobahrie A., Ibanez L., Bullitt E., Aylward S. "Vessel
%   Enhancing Diffusion Filter", Insight Journal,
%   2007. http://hdl.handle.net/1926/558.
%
%   B has the same size and class as A.
%
%   Note that A should have a signed type (e.g. int8, int16,
%   double...). Images with unsigned types (e.g. uint8) will cause
%   intermediate results that should be negative to be truncated to
%   0, and the result will not be meaningful.
%
%   SIGMAMIN, SIGMAMAX are scalars with the limits of the multiscale
%   scheme, in the same units as the image. They should be set to
%   roughly the diameters of the smallest and largest vessels in the
%   image. By default, SIGMAMIN=0.2, SIGMAMAX=2.0.
%
%   NUMSIGMASTEPS is a scalar with the number of scales for the
%   analysis. The scales change exponentially, not linearly. Casual
%   testing suggests that the final result does not depend heavily on
%   this parameter. By default, NUMSIGMASTEPS=10.
%
%   NUMITERATIONS is a scalar with the number of times the multiscale
%   anisotropic difussion method is run. In practice, a higher number
%   of iterations means more blurring along the vessels, which is
%   usually desirable. The result will depend heavily on the number
%   of iterations chosen. By default, NUMITERATIONS=10.
%
%   WSTRENGTH is a scalar that indicates the strength of anisotropic
%   diffusion. Casual testing suggests that the result doesn't depend
%   much on this value. By default, WSTRENGTH=25.0.
%
%   SENSITIVITY is a scalar that indicates the sensitivity to the
%   vesselness response. Casual testing suggests that the result
%   doesn't depend much on this value. By default, SENSITIVITY=5.0.
%
%   TIMESTEP is a scalar with the time step size in the difussion
%   process. It needs to be small enough to avoid divergence, but
%   otherwise casual testing suggests that the result doesn't depend
%   much on this value. For 3D images, TIMESTEP < 0.0625. By default,
%   TIMESTEP=0.001.
%
%   EPSILON is a scalar. It's a small number to ensure the positive
%   definiteness of the diffusion tensor. By default, EPSILON=0.01.
 
% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.4.5
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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

error('MEX file not found')
