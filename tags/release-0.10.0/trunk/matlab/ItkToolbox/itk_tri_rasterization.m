function varargout = itk_tri_rasterization(varargin)
% ITK_TRI_RASTERIZATION  Rasterization of triangular mesh to binary
% segmentation
%
% This function provides a Matlab interface to
% itk::TriangleMeshToBinaryImageFilter.
%
% BW = itk_tri_rasterization(TRI, X, RES, SIZE, ORIGIN)
%
%   TRI, X describe a triangular mesh.
%
%   TRI is a 3-column matrix. Each row represents the indices of the three
%   vertices that form a triangle. TRI as a whole represents the closed
%   surface.
%
%   X is a 3-column matrix. Each row represents the Cartesian coordinates
%   of a vertex on the surface, indexed by TRI values.
%
%   RES is a 3-vector with the voxel size in (row, column, slice) format.
%
%   SIZE is a 3-vector with the number of voxels in (row, column, slice)
%   format.
%
%   ORIGIN is a 3-vector with the coordinates of the centre of the
%   bottom-left image voxel, in (x, y, z) format.
%
%   BW is the output uint8 binary segmentation. Voxels inside the mesh will
%   be set to 1, and voxels outside to 0. Note that voxels with their
%   centre exactly on the mesh offer different results depending on whether
%   they are at the left, right, top or bottom of the mesh. This is a
%   limitation in itk::TriangleMeshToBinaryImageFilter. See
%   test_itk_tri_rasterization.m in Gerardus for an example.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
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
