function zi = mba_surface_interpolation(x, y, z, xi, yi)
% MBA_SURFACE_INTERPOLATION  Scattered data Multilevel B-spline interpolation
%
% This MEX-function uses the MBA library [1] to compute a Multilevel
% B-spline interpolated surface from a scattered set of points.
%
% ZI = MBA_SURFACE_INTERPOLATION(X, Y, Z, XI, YI)
%
%   X, Y, Z are column vectors of the same size with the 3D
%   coordinates of a scattered set of points.
%
%   XI, YI are column vectors of the same size (but can have a
%   different size from X, Y and Z) with the 2D coordinates of the
%   locations that we want to interpolate.
%
%   ZI is a vector of the same length as XI and YI, with the
%   interpolated values.
%
% [1] MBA - Multilevel B-Spline Approximation
% Library. http://www.sintef.no/Projectweb/Geometry-Toolkits/MBA/
%
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

% To compile this MEX-file in a 64-bit linux architecture, run from the
% gerardus/matlab directory
%
% >> mex -v -largeArrayDims -outdir PointsToolbox/ -f ./engopts.sh PointsToolbox/mba_surface_interpolation.cpp
%
% To compile in a 32-bit linux architecture, run (untested)
%
% >> mex -v -outdir PointsToolbox/ -f ./engopts.sh PointsToolbox/mba_surface_interpolation.cpp
%
% For Windows, Mac or other architectures, the corresponding section
% in the ../engopts.sh file will need to be edited before it
% compiles. I cannot test them, so I haven't touched those sections.
error('Compiled MEX function has not been found')
