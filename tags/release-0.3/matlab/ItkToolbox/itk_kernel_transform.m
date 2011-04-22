function itk_kernel_transform
% ITK_KERNEL_TRANSFORM  ITK warps between 3D point sets with a known
% correspondence
%
% This MEX-function provides a Matlab interface to run ITK's Kernel
% Transforms:
%
%   itk::ElasticBodySplineKernelTransform
%   itk::ElasticBodyReciprocalSplineKernelTransform
%   itk::ThinPlateSplineKernelTransform
%   itk::ThinPlateR2LogRSplineKernelTransform
%   itk::VolumeSplineKernelTransform
%
% YI = ITK_KERNEL_TRANSFORM(X, Y, XI, TYPE)
%
%   X, Y are 3-column matrices with N rows. Each row has the
%   coordinates of a point. The warp is defined so that
%   X(i,:)->Y(i,:).
%
%   XI is a 3-column matrix with M rows. Each row has the coordinates
%   of a point to be warped.
%
%   YI has the same dimensions as XI. YI contains the coordinates of
%   the warped points.
%
%   TYPE is a string that allows to select the type of warp (no defaults):
%
%   'elastic':  itk::ElasticBodySplineKernelTransform
%   'elasticr': itk::ElasticBodyReciprocalSplineKernelTransform
%   'tps':      itk::ThinPlateSplineKernelTransform
%   'tpsr2':    itk::ThinPlateR2LogRSplineKernelTransform
%   'volume':   itk::VolumeSplineKernelTransform
%
%  To compile and install this function to PointToolbox, you need to have
%  installed the Insight Toolkit C++ files and Matlab.
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
%
% If cmake throws an error because it cannot find Matlab, then edit
% gerardus/matlab/CMakeLists.txt, and where it says
%
%    SET(MATLAB_ROOT "/usr/local/matlab/R2010b/")
%
%  change to your own Matlab root path.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

error('MEX file not found')
