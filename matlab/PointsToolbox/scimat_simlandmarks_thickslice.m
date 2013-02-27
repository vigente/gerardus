function [x, scimat] = scimat_simlandmarks_thickslice(scimat, p, dim, s)
% SCIMAT_SIMLANDMARKS_THICKSLICE  Simulate a human placing landmarks on the
% boundary of a "thick slice" segmentation
%
% This function simulates a human placing landmarks on a segmentation when
% the slices are thick or when only one out of every n slices will be
% segmented. Each slice will have a guaranteed number of points randomly
% selected.
%
% [X, SCIMAT2] = scimat_simlandmarks_thickslice(SCIMAT, P, DIM, S)
%
%   SCIMAT is a struct with a binary segmentation.
%
%   P is a scalar with the proportion of points from the boundary of each
%   slice that will be randomly selected. By default, P=0.1 (10% of the
%   points).
%
%   DIM is the dimension along which the 2D boundaries are computed and
%   sampled. By default, DIM=3.
%
%   S is a vector with the indices of the slices that will be sampled.

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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% check arguments
narginchk(1, 4);
nargoutchk(0, 2);

% defaults
if (nargin < 2)
    p = .1;
end
if (nargin < 3)
    dim = 3;
end
if (nargin < 4)
    s = 1:size(scimat.data, dim);
end

% size of the volume
sz = size(scimat.data);

% the size of the slice depends on which dimension we are scrolling through
conn = zeros(3, 3, 3);
switch dim
    case 1
        conn(2, :, :) = 1;
        slice = zeros(sz(3), sz(2));
    case 2
        conn(:, 2, :) = 1;
        slice = zeros(sz(3), sz(1));
    case 3
        conn(:, :, 2) = 1;
        slice = zeros(sz(1), sz(2));
end

% to find a good perimeter, we need to go slice by slice. This also mimicks
% what a human operator would do
im = bwperim(scimat.data, conn);

% clear up the input segmentation so that we can place here the output
% sampling points
scimat.data = zeros(sz);

% in each slice, keep only 100*p% of the points, to simulate a human
% placing landmarks
for S = s
    
    % get slice
    switch dim
        case 1
            slice = im(S, :, :);
        case 2
            slice = im(:, S, :);
        case 3
            slice = im(:, :, S);
    end
    
    % indices of boundary points
    idx = find(slice);
    
    % number of points in this slice
    Ptot = length(idx);
    
    % number of points to keep in this slice
    P = round(Ptot * p);
    
    % sample points from the boundary
    idx = idx(randperm(Ptot, P));
    
    % clear up the slice and put the sample points back
    slice(:) = 0;
    slice(idx) = 1;

        % get slice
    switch dim
        case 1
            scimat.data(S, :, :) = slice;
        case 2
            scimat.data(:, S, :) = slice;
        case 3
            scimat.data(:, :, S) = slice;
    end
    
end

% get coordinates of landmarks
[r, c, s] = ind2sub(sz, find(scimat.data));
x = scinrrd_index2world([r, c, s], scimat.axis);
