function skeleton_plot(cc, lab)
% SKELETON_PLOT  Plot a labelled segmentation and its skeleton
%
% SKELETON_PLOT(CC, LAB)
%
%   CC and LAB are the outputs from function skeleton_label().
%
%   This function plots the branches from the skeleton with the two end
%   points of each one (note than several branches can share the same end
%   point). If LAB is provided and not empty, the branches are plotted on
%   the labelled segmentation.
%
% See also: skeleton_label.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.1.0
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
error(nargchk(1, 2, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% defaults
col = 'b';

% save hold status for later
HOLD = ishold;

% plot
if (nargin > 1 && ~isempty(lab))
    imagesc(lab)
    axis xy
    hold on
    col = 'w';
end

% loop every branch in the skeleton
for I = 1:cc.NumObjects
    % list of voxels in current branch
    br = cc.PixelIdxList{I};
    
    % plot skeleton branches
    if (length(cc.ImageSize) == 2)
        % convert image linear indices to image r, c, s indices
        [r, c] = ind2sub(cc.ImageSize, br);
        
        % plot branch in order
        plot(c, r, col)
        hold on
        
        % plot starting and end points
        plot(c([1 end]), r([1 end]), [col 'o'])
    elseif (length(cc.ImageSize) == 3)
        % convert image linear indices to image r, c, s indices
        [r, c, s] = ind2sub(cc.ImageSize, br);
        
        % plot branch in order
        plot3(c, r, s, 'w')
        hold on
        
        % plot starting and end points
        plot(c([1 end]), r([1 end]), s([1 end]), [col 'o'])
    else 
        error('Can only print 2D or 3D skeletons')
    end
end

% recover hold status
if HOLD
    hold on
else
    hold off
end
