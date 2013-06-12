% test_scimat_dmatrix_imgeodesic.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010-2013 University of Oxford
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Toy example, ellipsoidal shape in 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ellipse axes
a = 20;
b = 9;
phi=45;

% generate ellipse as a binary mask
[gx, gy] = meshgrid(-31:31, -21:21);
im = ((gx)/a).^2+((gy)/b).^2 <= 1;
im = uint8(im);

% increase the difference between dark and bright areas, otherwise, the
% bright area will not be necessarily avoided
im(im ~= 0) = 100;

% give some intensity value to the background. Otherwise, if the cost on
% the background is zero, any path can be chosen, not necessarily one that
% goes close to the edge
im = im + 10;

% convert to SCI MAT format
scimat = scinrrd_im2nrrd(im, [1 1 1], [.5 .5 .5]);

% plot ellipse binary mask
hold off
imagesc([.5 .5 .5], [.5 .5 .5], scimat.data)
axis equal xy

% choose 4 points on the edge of the ellipse
x = [
    31.5 31 0
    31.5 12 0
    11 21.5 0
    52 21.5 0
    ];

% real world coordinates => indices
idx = scinrrd_world2index(x, scimat.axis);

% round to nearest voxel centre
idx = round(idx);

% recompute coordinates of the points
x = scinrrd_index2world(idx, scimat.axis);

% plot points
hold on
plot(x(:, 1), x(:, 2), 'ow')

% compute livewire geodesic distance
[d2, idx, p, d] = scimat_dmatrix_imgeodesic(scimat, idx);

% loop every landmark
for I = 1:size(x, 1)
    
    % plot ellipse binary mask
    hold off
    imagesc([.5 .5 .5], [.5 .5 .5], scimat.data)
    axis equal xy
    hold on
    plot(x(:, 1), x(:, 2), 'ow')
    plot(x(I, 1), x(I, 2), 'og')

    for J = 1:size(x, 1)
        
        % skip path from node to itself
        if (I==J)
            continue
        end
        
        % recover the shortest paths from the current point (source) to the
        % target point
        pth = graphpred2path(p(I, :), idx(J));
        
        % convert linear indices to r, c, s
        [r, c, s] = ind2sub(size(scimat.data), pth);
        
        % obtain x, y, z coordinates
        pth = scinrrd_index2world([r' c' s'], scimat.axis);
        
        % plot path
        plot(pth(:, 1), pth(:, 2), 'w')
        
    end
    
    pause
    
end
