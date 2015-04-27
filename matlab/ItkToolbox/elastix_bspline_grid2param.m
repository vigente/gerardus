function v = elastix_bspline_grid2param(varargin)
% ELASTIX_BSPLINE_GRID2PARAM  Convert control points from grid to
% TransformParameters vector form.
%
% V = ELASTIX_BSPLINE_GRID2PARAM(GX, GY)
%
%   GX, GY are matrices describing a grid of B-spline control points
%   coordinates or displacements. This is a similar format to the output of
%   NDGRID. They follow the Matlab convention that X->columns and Y->rows.
%
%   Important: Function elastix_bspline_grid produces control point
%   coordinates, not displacements. So if you apply
%   ELASTIX_BSPLINE_GRID2PARAM to the output of elastix_bspline_grid, you
%   will get a vector of coordinates, not of displacements. Note that
%   elastix uses displacements in TransformParameters. To obtain the
%   displacements, see example below.
%
%   V is a vector with all the values in GX and GY formatted for an elastix
%   TransformParameters transform field.
%
% V = ELASTIX_BSPLINE_GRID2PARAM(GX, GY, ..., GN)
% 
%   Syntax for N-dimensional grids.
%
% Example:
%
%   % assume we have a B-spline transform t
%
%   % coordinates of displaced control points
%   [gx, gy] = elastix_bspline_grid(t);
%
%   % coordinates of control points without displacement
%   t0 = t;
%   t0.TransformParameters(:) = 0;
%   [gx0, gy0] = elastix_bspline_grid(t0);
%
%   % convert the displacements back from grid to vector format
%   v = elastix_bspline_grid2param(gx - gx0, gy - gy0);
%
%   % check that we have recovered the original transform parameters vector
%   % (an output of 0 means no error)
%   any(abs(v - t.TransformParameters) > 1e-12)
%
% See also: elastix_bspline_grid.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2015 University of Oxford
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% check input arguments
narginchk(1, Inf);
nargoutchk(0, 1);

% number of dimensions
D = length(varargin);

% swap rows<->columns from Matlab to elastix convention
varargin = cellfun(@(x) permute(x, [2 1 3:D]), varargin, ...
    'UniformOutput', false);

% linearize and concatenate in elastix format
varargin = cellfun(@(x) x(:)', varargin, ...
    'UniformOutput', false);
v = [varargin{:}];
