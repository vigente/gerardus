% ITK_PSTRANSFORM: Spatial transformation (i.e. warp) of a set of points,
% defined from a known landmark correspondence
%
% This MEX-function provides a Matlab interface to run ITK's Kernel
% and B-spline transforms:
%
%   itk::ElasticBodySplineKernelTransform
%   itk::ElasticBodyReciprocalSplineKernelTransform
%   itk::ThinPlateSplineKernelTransform
%   itk::ThinPlateR2LogRSplineKernelTransform
%   itk::VolumeSplineKernelTransform
%   itk::BSplineScatteredDataPointSetToImageFilter
%
%
% YI = ITK_PSTRANSFORM(TRANSFORM, X, Y, XI)
%
%   X, Y are 2-column (2D) or 3-column (3D) matrices. Each row has
%   the coordinates of a point. The warp is defined so that
%   X(i,:)->Y(i,:).
%
%   XI is a matrix with the same number of columns as X, Y. Each row
%   has the coordinates of a point to be warped.
%
%   YI has the same dimensions as XI. YI contains the coordinates of
%   the warped points.
%
%   TRANSFORM is a string that allows to select the type of warp (no
%   defaults):
%
% YI = ITK_PSTRANSFORM('elastic', X, Y, XI)
% YI = ITK_PSTRANSFORM('elasticr', X, Y, XI)
% YI = ITK_PSTRANSFORM('tps', X, Y, XI)
% YI = ITK_PSTRANSFORM('tpsr2', X, Y, XI)
% YI = ITK_PSTRANSFORM('volume', X, Y, XI)
%
%   'elastic':  itk::ElasticBodySplineKernelTransform
%   'elasticr': itk::ElasticBodyReciprocalSplineKernelTransform
%   'tps':      itk::ThinPlateSplineKernelTransform
%   'tpsr2':    itk::ThinPlateR2LogRSplineKernelTransform
%   'volume':   itk::VolumeSplineKernelTransform
%
% YI = ITK_PSTRANSFORM('bspline', X, Y, XI, ORDER, LEVELS)
%
%   'bspline':  itk::BSplineScatteredDataPointSetToImageFilter
%
%   By Nicholas J. Tustison, James C. Gee in the Insight Journal
%   paper: http://hdl.handle.net/1926/140
%
%   ORDER is an integer with the B-spline order. By default, ORDER=3,
%   and the B-spline is cubic.
%
%   LEVELS is an integer with the number of multi-resolution levels
%   in the algorithm. A higher number of levels will make the spline
%   more flexible and match the landmarks better. By default, LEVELS=5.
%
% This function must be compiled before it can be used from Matlab.
% If Gerardus' root directory is e.g. ~/gerardus, type from a
% linux shell
%
%    $ cd ~/gerardus/matlab
%    $ mkdir bin
%    $ cd bin
%    $ cmake ..
%    $ make install

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

error('MEX file not found')
