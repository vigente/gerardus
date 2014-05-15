function h = quadsurf(x, y, z, c, varargin)
% QUADSURF  Quadrangular mesh surface plot
%
% quadsurf(X, Y, Z)
%
%   X, Y, Z are matrices with the coordinates of the mesh vertices. Note
%   that with this function it is possible to plot meshes that bend on
%   themselves or are rolled up. This is not possible with Matlab's
%   function surf(), that requires a unique value of Z for each (X,Y) pair.
%
%   This is an example of matrix X for a mesh with 10 rows and 7 columns.
%
%   X(1,1)-----X(1,2)---...---X(1,7)
%    |          |        |      |
%    |          |        |      |
%    |          |        |      |
%   X(2,1)-----X(2,2)---...---X(2,7)
%    |          |        |      |
%   ...        ...      ...    ...
%    |          |        |      |
%   X(10,1)----X(20,2)--...---X(20,7)
%
% quadsurf(X, Y, Z, C, <parameter/value pairs>)
%
%   C is a matrix of the same size as X, Y and Z, and specifies the colour
%   of the vertices by indexing into the colormap. By default, C=Z.
%
%   The X,Y,Z,C quad can be followed by parameter/value pairs to specify
%   additional properties of the Patch, e.g. 
%
%     quadsurf(x, y, z, c, 'EdgeColor', 'white');
%
% H = quadsurf(...)
%
%   H returns a patch handle. Patches are children of AXES objects.
%
% See also: trisurf, trimesh, surf, mesh.

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

% check arguments
narginchk(3, Inf);
nargoutchk(0, 1);

% defaults
if (nargin < 4)
    c = z;
end

% size of the inputs
[R, C] = size(x);
if (any(size(y)~=[R, C]) || any(size(z)~=[R, C]))
    error('X, Y, Z must have the same size')
end

% patch takes quadrangles as rows of 4 elements. Here we create matrices
% with the same size as X,Y,Z, for each of the 4 vertices of a quadrangle.
% Note that we have a lot of redundancy, but this is the patch convention
%
% vx1       vx2
%  ----------
%  |        |
%  |        |
%  |        |
%  |        |
%  ----------
% vx4       vx3
vx1 = x(1:R-1, 1:C-1);
vx2 = x(1:R-1, 2:C);
vx3 = x(2:R, 2:C);
vx4 = x(2:R, 1:C-1);
vy1 = y(1:R-1, 1:C-1);
vy2 = y(1:R-1, 2:C);
vy3 = y(2:R, 2:C);
vy4 = y(2:R, 1:C-1);
vz1 = z(1:R-1, 1:C-1);
vz2 = z(1:R-1, 2:C);
vz3 = z(2:R, 2:C);
vz4 = z(2:R, 1:C-1);
c1 = c(1:R-1, 1:C-1);
c2 = c(1:R-1, 2:C);
c3 = c(2:R, 2:C);
c4 = c(2:R, 1:C-1);

% format the vertices matrices, and create the patch image
h = patch([vx1(:)'; vx2(:)'; vx3(:)'; vx4(:)'], ...
    [vy1(:)'; vy2(:)'; vy3(:)'; vy4(:)'], ...
    [vz1(:)'; vz2(:)'; vz3(:)'; vz4(:)'], ...
    [c1(:)'; c2(:)'; c3(:)'; c4(:)'], varargin{:});
