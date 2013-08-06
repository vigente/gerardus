function varargout = cgal_fixed_alpha_shape3(varargin)
% CGAL_FIXED_ALPHA_SHAPE3  Individual alpha-shapes of a 3D set of points
%
% TRI = cgal_fixed_alpha_shape3(X, ALPHA)
%
%   X is a 3-column matrix. X(i, :) contains the xyz-coordinates of a 3D
%   point.
%
%   ALPHA is a vector of scalar alpha values, alpha=R^2, where R is the
%   probe radius.
%
%   TRI is a cell array of the same length as ALPHA. Cell TRI{i} contains
%   the alpha shape triangulation for ALPHA{i}. Each row contains the 3
%   nodes that form one triangular facet in the mesh. The i-th mesh can be
%   visualised running:
%
%     >> trisurf(tri{i}, x)
%
%
% This function uses CGAL's implementation of fixed alpha shapes [1]. Fixed
% alpha shapes are more efficient when only the shape for one or a few
% alpha values is required. When many alpha values are required, it may be
% faster to use cgal_alpha_shape3().
%
% However, note that Matlab function alphavol() implemented by Jonas
% Lundgren and provided as a third-party function in Gerardus seems to be
% faster than either of the CGAL MEX functions, at least for a single alpha
% value.
%
% [1] http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Alpha_shapes_3/Chapter_main.html
%
% See also: alphavol, cgal_alpha_shape3, scimat_lconvhull_smoothing

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
