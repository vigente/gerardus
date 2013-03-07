function a = im2imat(im)
% IM2IMAT  Sparse distance matrix between 2D or 3D image voxels, weighted
% by voxel intensities
%
% A = IM2IMAT(IM)
%
%   IM is an image volume with dimensions (R, C, S).
%
%   A is a sparse matrix with dimensions (R*C*S, R*C*S), where element
%   (i,j) is the mean intensity between voxels with linear indices i and j.
%
%   Voxels with an Inf intensity are skipped.
%
% ... = IM2IMAT(..., RES) [This option only available in the MEX version]
%
%   RES is a row vector with the voxel size of [row, column, slice] (2D) or
%   [row, column, slice] (3D). By default, RES=[1.0 1.0 1.0].
%
%   This function has a slow Matlab implementation (using loops), but a
%   fast MEX version is provided with Gerardus too.
%
% See also: seg2dmat.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010-2013 University of Oxford
% Version: 0.1.2
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

warning('Warning: Running Matlab version, slower than compiled MEX version')

a = sparse([],[],[], numel(im), numel(im), 26*length(find(im)));
for S = 1:size(im, 3),
    for C = 1:size(im, 2),
        for R = 1:size(im, 1)
            
            % linear index of current voxel
            idx = sub2ind(size(im), R, C, S);

            % if current voxels is Inf, we don't need to include it in the
            % graph explicitly
            if (~isinf(im(idx)))
                
                % indices of the 27 voxels forming a cube around the
                % current voxel, incluiding itself
                [rg, cg, sg] = ndgrid(...
                    max(1, R-1):min(size(im,1), R+1), ...
                    max(1, C-1):min(size(im,2), C+1), ...
                    max(1, S-1):min(size(im,3), S+1));
                
                % linear indices of the 27 cube voxels
                cube = sub2ind(size(im), rg(:), cg(:), sg(:))';
                
                % don't connect the central voxel to itself
                cube = cube(cube ~= idx);
                
                % indices of cube voxels that are not Inf
                cube = cube(~isinf(im(cube)));
                
                % the weight of the edge between the current voxel and each
                % connected voxel is the mean intensity between both voxels
                a(idx, cube) = (im(idx) + im(cube))/2;
                a(cube, idx) = (im(idx) + im(cube))/2;
            
            end % if (~isinf(im(idx)))
            
        end
    end
end
