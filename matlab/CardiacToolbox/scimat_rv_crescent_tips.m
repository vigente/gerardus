function [tip1, tip2] = scimat_rv_crescent_tips(scimat, m)
% SCIMAT_RV_CRESCENT_TIPS  Extract the tips of the crescent-shaped curve
% in all slices of the Right Ventricle.
%
% [X1, X2] = scimat_rv_crescent_tips(SCIMAT, M)
%
%   X1, X2 are 3-colum matrices where each row are the real world
%   coordinates of each of the tips of the RV crescent shape.
%
%   SCIMAT is the struct with the RV segmentation mask (see "help scimat"
%   for details).
%
%   M is a 3-column matrix, where each row has the real world coordinates
%   of a centroid (typically, the central Left Ventricle's central curve
%   computed with scimat_centroids()). There must be as many centroids as
%   slices. Centroids with NaN coordinates will be skipped.
%
%   In the case of a RV slice that has no corresponding LV centroid, the
%   closest LV centroid will be used.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010,2014 University of Oxford
% Version: 0.2.0
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
narginchk(2, 2);
nargoutchk(0, 2);

% volume size
sz = [scimat.axis.size];

% if (size(m, 1) ~= sz(3))
%     error('There must be a centroid per slice in SCIMAT, even if the centroid in NaN')
% end
if (size(m, 2) ~= 3)
    error('M must be a 3-column matrix')
end

% init outputs
tip1 = nan(sz(3), 3);
tip2 = nan(sz(3), 3);

% iterate slices
for I = 1:sz(3)
    % extract slice
    im = scimat.data(:,:, I);
    
    % get linear indices of all pixels in current slice
    idx = find(im);
    
    % if RV slice is empty, skip it
    if isempty(idx)
        continue
    end
    
    % convert linear index to multiple subscripts
    [ir, ic] = ind2sub( sz(1:2), idx );
    
    % convert indices to real world coordinates and make colum vectors
    x = scimat_index2world( [ ir, ic, I+zeros(length(ir), 1) ], ...
        scimat.axis );
    
    % compute a centroid for the slice
    xm = mean(x, 1);
    
    % compute the distance from the slice centroid to each axis centroid
    d = dmatrix(xm', m');
    
    % find closest axis centroid to the slice centroid
    [~, inn] = min(d);
    
    % center the slice points around the axis centroid
    x = x - m(inn * ones(size(x, 1), 1), :);
    
    % compute polar coordinates of the centered pixel coordinates
    [phi, th] = cart2pol(x(:,1), x(:,2));
    
    % assume that the point with the largest and smallest azimuth values
    % are the crescent tips
    [~, idx1] = max(phi);
    [~, idx2] = min(phi);
    
    % assign the tip points to the output variable, undoing the previous
    % centering
    tip1(I, :) = x(idx1, :) + m(inn, :);
    tip2(I, :) = x(idx2, :) + m(inn, :);
        
end
