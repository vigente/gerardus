function [tip1, tip2] = scinrrd_rv_crescent_tips(nrrd, m)
% SCINRRD_RV_CRESCENT_TIPS  Extract the tips of the crescent-shaped curve
% in all slices of the Right Ventricle
%
% [X1, X2] = SCINRRD_RV_CRESCENT_TIPS(NRRD, M)
%
%   X1, X2 are 3-colum matrices where each row are the real world
%   coordinates of each of the tips of the RV crescent shape.
%
%   NRRD is the struct with the RV segmentation mask.
%
%   M is a 3-column matrix, where each row has the real world coordinates
%   of a centroid (typically, the central Left Ventricle's central curve
%   computed with scinrrd_centroids()). There must be as many centroids as
%   slices. Centroids with NaN coordinates will be skipped.
%
%   In the case of a RV slice that has no corresponding LV centroid, the
%   closest LV centroid will be used.
%
%
%   Note on SCI NRRD: Software applications developed at the University of
%   Utah Scientific Computing and Imaging (SCI) Institute, e.g. Seg3D,
%   internally use NRRD volumes to store medical data.
%
%   When label volumes (segmentation masks) are saved to a Matlab file
%   (.mat), they use a struct called "scirunnrrd" to store all the NRRD
%   information:
%
%   >>  scirunnrrd
%
%   scirunnrrd = 
%
%          data: [4-D uint8]
%          axis: [4x1 struct]
%      property: []

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010 University of Oxford
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
error( nargchk( 2, 2, nargin, 'struct' ) );
error( nargoutchk( 0, 2, nargout, 'struct' ) );

% volume size
sz = [nrrd.axis.size];

% if (size(m, 1) ~= sz(3))
%     error('There must be a centroid per slice in NRRD, even if the centroid in NaN')
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
    im = nrrd.data(:,:, I);
    
    % get linear indices of all pixels in current slice
    idx = find(im);
    
    % if RV slice is empty, skip it
    if isempty(idx)
        continue
    end
    
    % convert linear index to multiple subscripts
    [ir, ic] = ind2sub( sz(1:2), idx );
    
    % convert indices to real world coordinates and make colum vectors
    x = scinrrd_index2world( [ ir, ic, I+zeros(length(ir), 1) ], ...
        nrrd.axis );
    
    % compute a centroid for the slice
    xm = mean(x, 1);
    
    % compute the distance from the slice centroid to each axis centroid
    d = dmatrix(xm', m');
    
    % find closest axis centroid to the slice centroid
    [foo, inn] = min(d);
    
    % center the slice points around the axis centroid
    x = x - m(inn * ones(size(x, 1), 1), :);
    
    % compute polar coordinates of the centered pixel coordinates
    [phi, th] = cart2pol(x(:,1), x(:,2));
    
    % assume that the point with the largest and smallest azimuth values
    % are the crescent tips
    [foo, idx1] = max(phi);
    [foo, idx2] = min(phi);
    
    % assign the tip points to the output variable, undoing the previous
    % centering
    tip1(I, :) = x(idx1, :) + m(inn, :);
    tip2(I, :) = x(idx2, :) + m(inn, :);
        
end
