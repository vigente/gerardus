function im = itk_imfilter(~, ~)
% ITK_IMFILTER: Run ITK filter on a 2D, 3D, 4D or 5D image
%
% This MEX function is a multiple-purpose wrapper to be able to run
% all ITK filters that inherit from itk::ImageToImageFilter on a
% Matlab 2D image or 3D, 4D or 5D image volume.
%
% B = ITK_IMFILTER(TYPE, A, [FILTER PARAMETERS])
%
%   TYPE is a string with the filter we want to run. See below for a whole
%   list of options.
%
%   A is a 2D matrix or 3D, 4D or 5D volume with the image or
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
%   A can also be a SCI MAT struct, A = scimat, with the following fields:
%
%     scimat.data: 2D or 3D array with the image or segmentation, as above
%     scimat.axis: 3x1 struct array with fields:
%       scimat.axis.size:    number of voxels in the image
%       scimat.axis.spacing: voxel size, image resolution
%       scimat.axis.min:     real world coordinates of image origin
%       scimat.axis.max:     ignored
%       scimat.axis.center:  ignored
%       scimat.axis.label:   ignored
%       scimat.axis.unit:    ignored
%
%   (An SCI MAT struct is the output of Matlab's function scinrrd_load(),
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
% [B, V, W] = itk_imfilter('dandist', A).
%
%   (itk::DanielssonDistanceMapImageFilter) Compute unsigned distance map
%   for a binary mask. Distance values are given in voxel coordinates.
%
%   B has the same size as A and type double. Each element in B
%   contains an approximation to the Euclidean distance of that voxel
%   to the closest foreground voxel, in index units.
%
%   V has the same size and type as A. V is a Voronoi partition of A,
%   using the same indices.
%
%   W has size (3,R,C,S) if A has size (R,C,S), and type int64. Each
%   3-vector W(:,i,j,k) is a vector pointing to the closest
%   foreground voxel from A(i,j,k).
%
% B = ITK_IMFILTER('maudist', A)
%
%   (itk::SignedMaurerDistanceMapImageFilter) Compute signed distance map
%   for a binary mask. Distance values are given in real world coordinates,
%   if the input image is given as a SCI MAT struct, or in voxel units, if
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
%   Note: A should have a signed type (e.g. int16, single). Images
%   with unsigned types (e.g. uint16) will cause intermediate results
%   that should be negative to be truncated to 0, and the result will
%   not be meaningful. The best compromise between accuracy and
%   saving memory seems to be type single (= float).
%
%   Note: While it is possible to run the filter on a SCI MAT struct,
%   results seem better if run directly on the image. The
%   filter doesn't seem to be spacing invariant.
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
%   ISSIGMASTEPLOG is a boolean that determines whether the
%   intermediate scales between SIGMAMIN to SIGMAMAX are distributed
%   logarithmically (true) or linearly (false). The latter seems to
%   work better for small ranges.
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
%
% B = itk_imfilter('hesves', A, SIGMAMIN, SIGMAMAX, NUMSIGMASTEPS, ISSIGMASTEPLOG)
%
%   (itk::MultiScaleHessianSmoothed3DToVesselnessMeasureImageFilter)
%   Vesselness measure from a multiscale scheme based on
%   eigenanalysis of the Hessian.
%
%   Enquobahrie A., Ibanez L., Bullitt E., Aylward S. "Vessel
%   Enhancing Diffusion Filter", Insight Journal,
%   2007. http://hdl.handle.net/1926/558.
%
%   B has the same size as A, but is always of type double.
%
%   Input arguments are the same as the four first input arguments of
%   filter "advess" above.
%
% B = itk_imfilter('median', A, RADIUS)
%
%   (itk::MedianImageFilter)
%   Median of a rectangular neighbourhood.
%
%   B has the same size and class as A.
%
%   RADIUS is a vector of scalars with the half-size of the filter's
%   box in each dimension. E.g. RADIUS=[2, 3, 4] means that the
%   median is computed in a rectangular neighbourhood of [5, 7, 9]
%   voxels.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.5.0
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
