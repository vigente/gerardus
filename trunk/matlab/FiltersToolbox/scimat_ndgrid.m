function [xg, yg, zg, tg] = scimat_ndgrid(scimat, ri, ci, si, fi)
% SCIMAT_NDGRID  Generation of arrays for 2D to 4D SCIMAT image volumes.
%
% [XG, YG, ZG] = scimat_ndgrid(SCIMAT,RI,CI,SI,FI)
%
%   SCIMAT is a struct with an image volume (see "help scimat" for
%   details). 
%
%   XG, YG, ZG are arrays with the x-, y-, z-coordinates of the voxels in
%   SCIMAT.
%
%   Note that XG values change with columns, and YG values change with
%   rows, to accommodate the usual coordinate frame convention.
%
% [XG, YG, ZG] = scimat_ndgrid(SCIMAT, RI, CI, SI, FI)
%
%   RI, CI, SI, FI are vectors of voxel indices (rows, columns, slices and frames). 
%   The output grid will be generated only for the corresponding image block.
%   For example, RI=1:4, CI=3:6, SI=5:7 means that the grid of coordinates
%   will be generated only for the block of voxels (1:4, 3:6, 5:7).
%
% See also: ndgrid, scimat_index2world.

% Author(s): Ramon Casero <rcasero@gmail.com>,
% Ben Villard <b.016434@gmail.com>
% Copyright © 2010-2015 University of Oxford
% Version: 0.4.1
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

% check arguments
D = length(scimat.axis);
narginchk(1, 5);
nargoutchk(0, D);

% simplify notation
NR = size(scimat.data, 1);
NC = size(scimat.data, 2);
NS = size(scimat.data, 3);
NF = size(scimat.data, 4);

% defaults 
if (nargin < 2 || isempty(ri))
    ri = 1:NR; % Creates indices of data length (rows).
end
if (nargin < 3 || isempty(ci))
    ci = 1:NC; % Creates indices of data length (columns).
end
if (nargin < 4 || isempty(si))
    si = 1:NS; % Creates indices of data length (number of slices).
end
if (nargin < 5 || isempty(fi))
    fi = 1:NF; % Creates indices of data length (number of time frames).
end

% generate 3D grid of coordinates. X and Y have to be swapped so that the
% grid has the right number of rows and columns
if (D == 2)
    
    [rg, cg] = ndgrid(ri, ci);
    xyzt = scimat_index2world([rg(:), cg(:)], scimat); % Converts indices to world coordinates
    xg = reshape(xyzt(:, 1), size(rg));
    yg = reshape(xyzt(:, 2), size(cg));
  
elseif (D == 3)
    
    [rg, cg, sg] = ndgrid(ri, ci, si);
    xyzt = scimat_index2world([rg(:), cg(:), sg(:)], scimat);% Converts indices to world coordinates
    xg = reshape(xyzt(:, 1), size(rg));
    yg = reshape(xyzt(:, 2), size(cg));
    zg = reshape(xyzt(:, 3), size(sg));
  
elseif (D == 4)
    [rg, cg, sg, fg] = ndgrid(ri, ci, si, fi);
    xyzt = scimat_index2world([rg(:), cg(:), sg(:), fg(:)], scimat);% Converts indices to world coordinates
    xg = reshape(xyzt(:, 1), size(rg));
    yg = reshape(xyzt(:, 2), size(cg));
    zg = reshape(xyzt(:, 3), size(sg));
    tg = reshape(xyzt(:, 4), size(fg));
    
else
    error('Wrong number of dimensions in scimat.data')
end


end
