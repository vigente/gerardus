% test_scimat_dmatrix_imgeodesic.m

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013-2014 University of Oxford
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

% plot ellipse
subplot(2, 1, 1)
hold off
imagesc([.5 62.5], [.5 42.5], im)
axis equal

% parameters for edge detector filter
sigma = [2 2];
maxerr = 0.01*ones(1, 2);

% run Canny edge filter
[~, imed] = itk_imfilter('canny', single(im), sigma.^2, [], [], maxerr);

% plot edge detection before thresholds
subplot(2, 1, 2)
imagesc(imed)

% apply thresholds
uppthr = .1;
lowthr = uppthr/2;
[imed, imed2] = itk_imfilter('canny', single(im), sigma.^2, uppthr, lowthr, maxerr);

% compute distance map (Maurer distance is much faster than Danielsson),
% and make all values positive
imd = abs(itk_imfilter('maudist', imed));

% plot distance map
imagesc(imd)

% convert to SCI MAT format
scimat = scinrrd_im2nrrd(imd, [1 1 1], [.5 .5 .5]);

% plot ellipse distance map binary mask
subplot(2, 1, 2)
hold off
imagesc([.5 62.5], [.5 42.5], scimat.data)
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
plot(x(:, 1), x(:, 2), 'or')

% compute livewire geodesic distance
[d2, idx, p, d] = scimat_dmatrix_imgeodesic(scimat, idx);

% loop every landmark
for I = 1:size(x, 1)
    
    % plot ellipse binary mask
    hold off
    subplot(2, 1, 2)
    imagesc([.5 62.5], [.5 42.5], scimat.data)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3D example, rat left ventricle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load downsampled MRI image or rat heart: Note this doesn't work because
% the papillary muscles create a low-cost shortcut between the septum and
% LV free wall
scimat = scimat_load('data/008-mri-downsampled-4.mha');

% instead, load a rough segmentation mask of the LV, with the papillary
% muscles removed with an alpha-shape
scimat = scimat_load('data/008-lvhull-downsampled-4.mha');

% plot one of the slices
subplot(2, 1, 1)
hold off
imagesc(scimat.data(:, :, 112))
colormap(gray)
axis equal

% to avoid numerical errors in the edge detector, normalize the voxel size
inc = min([scimat.axis.spacing]);
scimat.axis(1).spacing = scimat.axis(1).spacing / inc;
scimat.axis(2).spacing = scimat.axis(2).spacing / inc;
scimat.axis(3).spacing = scimat.axis(3).spacing / inc;

% parameters for edge detector filter
sigma = [2 2 2];

% CannyEdgeDetectionImageFilter only accepts input images with floating
% type (double or single)
scimat.data = single(scimat.data);

% run Canny edge filter
[~, imed2] = itk_imfilter('canny', scimat, sigma.^2);

% plot edge detection before thresholds
subplot(2, 1, 2)
hold off
imagesc(imed2(:, :, 112))
colormap(gray)
axis equal

% apply thresholds
uppthr = .1;
lowthr = uppthr/2;
imed = itk_imfilter('canny', scimat, sigma.^2, uppthr, lowthr);

% plot edge detection after thresholds
subplot(2, 1, 2)
hold off
imagesc(imed(:, :, 112))
colormap(gray)
axis equal

% compute distance map (Maurer distance is much faster than Danielsson).
% Make sure to make all distance values positive, otherwise Dijkstra will
% crash
scimat2.data = imed;
scimat2.axis = scimat.axis;
imd = abs(itk_imfilter('maudist', scimat2));

% plot distance map
subplot(2, 1, 2)
hold off
imagesc(imd(:, :, 112))
colormap(jet)
axis equal

% load landmarks on the LV endocardium
aux = scimat_load('data/008-lvhull-points.mha');
% aux = scimat_load('data/008-rvhull-points.mha');
[r, c, s] = ind2sub(size(aux.data), find(aux.data));
idx = [r c s];

% apply the same axis normalization we used above
aux.axis = scimat.axis;

% convert to real world coordinates
x = scinrrd_index2world(idx, scimat.axis);

% plot points
subplot(1, 1, 1)
hold off
plot3(x(:, 1), x(:, 2), x(:, 3), '.')
axis equal

% put the distance map in a SCI MAT struct, so that Euclidean distances
% between voxel centers are computed correctly
%
% increase distance values to increase contrast and penalise wandering away
% from an edge
scimat2.data = imd * 10;

% compute livewire geodesic distance
tic
[d, idx2, p] = scimat_dmatrix_imgeodesic(scimat2, idx);
toc

% create a segmentation mask to keep the paths
scimat3.data = zeros(size(scimat.data), 'uint8');
scimat3.axis = scimat.axis;

% loop every landmark
for I = 1:size(x, 1)
    
    % plot points
    subplot(1, 1, 1)
    hold off
    plot3(x(:, 1), x(:, 2), x(:, 3), '.')
    axis equal

    for J = 1:size(x, 1)
        
        % skip path from node to itself
        if (I==J)
            continue
        end
        
        % recover the shortest paths from the current point (source) to the
        % target point
        pth = graphpred2path(p(I, :), idx2(J));
        
        % record this path in the segmentation mask
        scimat3.data(pth) = 1;
        
        % convert linear indices to r, c, s
        [r, c, s] = ind2sub(size(scimat.data), pth);
        
        % obtain x, y, z coordinates
        pth = scinrrd_index2world([r' c' s'], scimat.axis);
        
        % plot path
        hold on
        plot3(pth(:, 1), pth(:, 2), pth(:, 3), 'r')
        
    end
    
    pause
    
end
