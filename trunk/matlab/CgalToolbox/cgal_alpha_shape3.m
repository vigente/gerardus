function varargout = cgal_alpha_shape3(varargin)
% CGAL_ALPHA_SHAPE3  Whole alpha-shape of a 3D set of points
%
% [ALPHALIM, TRI] = cgal_alpha_shape3(X, [], NCOMP)
%
%   X is a 3-column matrix. X(i, :) contains the xyz-coordinates of a 3D
%   point.
%
%   NCOMP is a scalar with the number of connected components in the output
%   alpha shape. By default, NCOMP=1.
%
%   ALPHALIM is a 3 vector:
%
%     ALPHALIM(1): minimum in the range of possible alpha values.
%     ALPHALIM(2): optimal alpha value for NCOMP connected  components.
%     ALPHALIM(3): maximum in the range of possible alpha values.
%
%   TRI is a cell containing a 3-column matrix with the alpha-shape
%   triangulation produced by ALPHALIM(2). This alpha-shape has NCOMP
%   connected components. Each row contains the 3 nodes that form one
%   triangular facet in the mesh. The mesh can be visualised running:
%
%     >> trisurf(tri{1}, x)
%
%
% [~, TRI] = cgal_alpha_shape3(X, ALPHA)
%
%   ALPHA is a vector of scalar alpha values, alpha=R^2, where R is the
%   probe radius. TRI is a cell array of the same length as ALPHA. Cell
%   TRI{i} contains the alpha shape triangulation for ALPHA{i}.
%
% This function uses CGAL's implementation of non-fixed alpha shapes [1].
% With non-fixed alpha shapes, the result for all alpha values is computed
% internally, and then only those the solutions requested by the user are
% extracted. If only an alpha value or a few are required, it may be faster
% to use cgal_fixed_alpha_shape3().
%
% However, note that Matlab function alphavol() implemented by Jonas
% Lundgren and provided as a third-party function in Gerardus seems to be
% faster than either of the CGAL MEX functions, at least for a single alpha
% value.
%
% [1] http://www.cgal.org/Manual/latest/doc_html/cgal_manual/Alpha_shapes_3/Chapter_main.html
%
% See also: alphavol, cgal_fixed_alpha_shape3, scimat_lconvhull_smoothing

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
% Version: 0.1.1
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
