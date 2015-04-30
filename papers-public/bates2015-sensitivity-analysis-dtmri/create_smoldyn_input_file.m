function [ ] = create_smoldyn_input_file( CSA, length, aspect_ratio, volume_fraction, diff_coeff, angle_range_per_voxel )
%CREATE_SMOLDYN_INPUT_FILE.M creates the text file which is the input to
%the Smoldyn software.
%
% Smoldyn: www.smoldyn.org/index.html
%
% This creates layers of cells at gradually increasing angles to simulate a
% section of myocardium. Cells are modelled as cuboids. Cells in each layer
% are parallel to each other.
%
% inputs:
%   CSA = cross-sectional area of cell, in um^2 (possible values for this study: 100, 167, 233, 300)
%   length = length of cell, in um (90, 110, 130, 150)
%   aspect_ratio = ratio of depth to width of cell cross-section (0.33, 0.56, 0.78, 1.00)
%   volume_fraction =  fraction of volume made of myocytes (0.75, 0.80, 0.85, 0.90)
%   diff_coeff = diffusivity of water in the myocardium, in um^2/ms (0.5, 1.17, 1.83, 2.50)
%   angle_range_per_voxel = (degrees) (30, 40, 50, 60)
%
% outputs:
%   output.txt file for use with the Smoldyn software
%
% notes:
%   positions data will be written to a folder called results, and are
%   called trajs1.txt ... trajs4.txt

% Author: Joanne Bates <joanne.bates@eng.ox.ac.uk>
% Copyright (c) 2015 University of Oxford
% Version: 0.1.0
% Date: 28 April 2015
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

              
% PARAMETERS
max_theta_voxel = 45 + (angle_range_per_voxel/2);
min_theta_voxel = 45 - (angle_range_per_voxel/2);

voxel_size = 600;       % um (incl the 200um boundary in each direction)
max_mols = 4000000;
no_mols = max_mols/4;   % divided into smaller parts, as can't make files too big
step_time = 0.1;        % ms
total_time = 10;        % ms

% calculated cell parameters
width = sqrt(CSA/aspect_ratio);
depth = CSA/width;

% calculation of gap between cells
% solve equation VF = (l*w*d)/((l+g)*(w+g)*(d+g))
% rearranged to cubic equation
% g^3 + (w+l+d)*g^2 + (wd + ld + wl)*g + (wdl - wdl/VF) = 0
% use roots to solve ax^3 + bx^2 + cx + d = 0
a = 1;
b = width + length + depth;
c = width*depth + width*length + length*depth;
d = (width*depth*length) - ((width*depth*length)/volume_fraction);
p = [a b c d];
r = roots(p);
r = r(imag(r)==0);
r = sort(r, 'descend');
gap = r(1);

% BASIC SET-UP
% geometry
distance_between_planes = depth + gap;
no_of_planes = floor(voxel_size/distance_between_planes);

% create the initial cornerpoint & the cell angles for each plane.
% initial cornerpoint randomly between 1 and width of cell to allow
% variable overlapping of cells from one plane to next.
starting_cornerpoints = zeros(no_of_planes, 2);
cell_angles = zeros(no_of_planes, 1);
for i = 1:(no_of_planes)
    starting_cornerpoints(i,:) = [rand(1,1)*width, (i-1)*distance_between_planes];
    cell_angles(i,:) = min_theta_voxel + (i-1)*((max_theta_voxel-min_theta_voxel)/(no_of_planes-1));
end

total_cells = ceil(no_of_planes * voxel_size * voxel_size / (length * width)); % just a guess to initialise the matrices
surfaces = zeros(108,total_cells);
cells = zeros(3,total_cells);
start_pt = 1;

clear a b c d i p r

% CALCULATIONS IN EACH PLANE
for k = 1:no_of_planes
    theta = cell_angles(k);
    z_position1 = starting_cornerpoints(k,2);   % edge of the plane
    z_position2 = z_position1 + depth;          % other edge of the plane
    start_corner(1,:) = [starting_cornerpoints(k,1),0];
    index = 1;
        
    no_of_columns = ceil(voxel_size/((width+gap)*max(cosd(theta), 0.1)));
    no_of_rows = ceil(voxel_size/((length+gap)*max(sind(theta), 0.1)));
    no_of_cells = no_of_columns*no_of_rows;
    corners = zeros(no_of_cells,2,4);
    % new start corner = old start corner + distance to top RH corner +
    % function of gap between cells
    for i = 2:no_of_rows
        start_corner(i,:) =  start_corner(i-1,:) + [-length*cosd(theta), length*sind(theta)] + [-gap*cosd(theta), gap*sind(theta)];
    end
    start_corner(no_of_rows+1,:) = start_corner(1,:) - [-length*cosd(theta), length*sind(theta)] - [-gap*cosd(theta), gap*sind(theta)];
    for i = no_of_rows+2:no_of_rows*2
        start_corner(i,:) =  start_corner(i-1,:) - [-length*cosd(theta), length*sind(theta)] - [-gap*cosd(theta), gap*sind(theta)];
    end
    for j = 1:no_of_rows*2
        a = start_corner(j,:);
        for m = 1:no_of_columns
            b = a + [width*sind(theta), width*cosd(theta)];     % top RH corner
            d = a + [-length*cosd(theta), length*sind(theta)];  % bottom LH corner
            c = d + [width*sind(theta), width*cosd(theta)];     % bottom RH corner
            corners(index,:,:) = [a;b;c;d]';
            a = b + [gap*sind(theta), gap*cosd(theta)];     % top LH corner of next cell in x direction
            index = index + 1;
         end
     end

    % for all cells where any corner is outside of the voxel, turn corners into
    % NaNs and then remove these so only have cells left which are fully in the
    % voxel.
    corners = reshape(corners, [no_of_cells*2, 8]);

    for i = 1:no_of_cells*2
        if max(corners(i,:)) >= voxel_size
            corners(i,:) = [nan nan nan nan nan nan nan nan];
        end
        if min(corners(i,:)) < 0
            corners(i,:) = [nan nan nan nan nan nan nan nan];
        end
    end

    corners = corners(isfinite(corners(:,1)),:);
    no_of_cells = size(corners, 1);
    corners = reshape(corners, [no_of_cells, 2, 4]);

    a = squeeze(corners(:,:,1));
    b = squeeze(corners(:,:,2));
    c = squeeze(corners(:,:,3));
    d = squeeze(corners(:,:,4));
    
    % create the required vectors for the Smoldyn text
    % z direction co-ordinates for the layer
    z_vector_1 = z_position1*ones(no_of_cells,1);
    z_vector_2 = z_position2*ones(no_of_cells,1);
    % each of the 8 corners of the cuboid. 1-4 are closest to z axis.
    corner1 = [a, z_vector_1];
    corner2 = [b, z_vector_1];
    corner3 = [c, z_vector_1];
    corner4 = [d, z_vector_1];
    corner5 = [a, z_vector_2];
    corner6 = [b, z_vector_2];
    corner7 = [c, z_vector_2];
    corner8 = [d, z_vector_2];
    % cuboid is made up of 12 triangles, 2 for each surface
    tri1 = [corner1 corner2 corner3];
    tri2 = [corner3 corner4 corner1];
    tri3 = [corner2 corner6 corner7];
    tri4 = [corner7 corner3 corner2];
    tri5 = [corner6 corner5 corner8];
    tri6 = [corner8 corner7 corner6];
    tri7 = [corner1 corner4 corner8];
    tri8 = [corner8 corner5 corner1];
    tri9 = [corner4 corner3 corner7];
    tri10 = [corner7 corner8 corner4];
    tri11 = [corner2 corner1 corner5];
    tri12 = [corner5 corner6 corner2];
    voxels_size = (voxel_size * ones(15,1))';
    voxels_size(2) = voxel_size;
    voxels_size(5) = voxel_size;
    voxels_size(7) = voxel_size;
    voxels_size(10) = voxel_size;
    voxels_size(13) = voxel_size;
    % midpoint of cuboid is mean of opposite corners
    cells_plane = (corner1 + corner7)/2;
    cells_plane = cells_plane';

    % create the matrices for Smoldyn text
    surfaces_plane = [tri1 tri2 tri3 tri4 tri5 tri6 tri7 tri8 tri9 tri10 tri11 tri12];
    surfaces_plane = surfaces_plane'; 

    % add the data for that plane into the full matrices
    total_cells = start_pt + no_of_cells - 1;
    surfaces(:,start_pt:total_cells) = surfaces_plane;
    cells(:,start_pt:total_cells) = cells_plane;
   
    start_pt = start_pt + no_of_cells;

end

clear axis_end axis_ends extra_in_x gap i j m n p q no_in_y no_in_x theta
clear midpoint_gap_x midpoint_gap_x_1 midpoint_gap_x_2 midpoint_gap_y_x midpoint_gap_y_y 
clear midpoint_top midpoint_top_1 midpoint_top_2 midpoints
clear bottom_end index mols no_of_cylinders no_points normal_vector radii top_end z_position z_vector

% full matrices for writing to the text file
index2 = (1:total_cells);
surfaces = surfaces(:,1:total_cells);
cells = cells(:,1:total_cells);
surfaces = [index2; surfaces];
cells = [index2; index2; cells];

% WRITE THE DATA TO THE SMOLDYN TEXT FILE
% open the basic file
fid = fopen('output.txt', 'a');  % appends the existing text

% write to the file
words1 = 'graphics opengl\ngraphic_iter 1\n\ndim 3\n\n#CSA %d\n#length %d\n#aspect ratio %d\n#volume fraction %d\n#diff coeff %d\n#angles per voxel %d\n\n';
fprintf(fid, words1, CSA, length, aspect_ratio, volume_fraction, diff_coeff, angle_range_per_voxel);
words2 = '\nboundaries 0 0 %d \nboundaries 1 0 %d \nboundaries 2 0 %d \nspecies one two three four\n\n';
fprintf(fid, words2, voxel_size, voxel_size, voxel_size);
words3 = '\nmax_mol %d \ndifc all %1.1f\ncolor all 1 0 0 \ndisplay_size all 1\n\n';
fprintf(fid, words3, max_mols, diff_coeff);
words4 = 'time_start 0 \ntime_stop %d \ntime_step %1.1f \n\n';
fprintf(fid, words4, total_time, step_time);
words5 = 'start_surface walls \naction both all reflect \ncolor both 0 0 0 \npolygon both edge \npanel rect +0 0 0 0 %d %d \npanel rect -0 %d 0 0 %d %d \npanel rect +1 0 0 0 %d %d \npanel rect -1 0 %d 0 %d %d \npanel rect +2 0 0 0 %d %d \npanel rect -2 0 0 %d %d %d \nend_surface\n\n';
fprintf(fid, words5, voxels_size);
words6 = 'start_surface surf%d \naction all both reflect \ncolor both purple 0.5 \nthickness 2 \npolygon front face \npolygon back edge \npanel tri %d %d %d %d %d %d %d %d %d \npanel tri %d %d %d %d %d %d %d %d %d \npanel tri %d %d %d %d %d %d %d %d %d \npanel tri %d %d %d %d %d %d %d %d %d \npanel tri %d %d %d %d %d %d %d %d %d \npanel tri %d %d %d %d %d %d %d %d %d \npanel tri %d %d %d %d %d %d %d %d %d \npanel tri %d %d %d %d %d %d %d %d %d \npanel tri %d %d %d %d %d %d %d %d %d \npanel tri %d %d %d %d %d %d %d %d %d \npanel tri %d %d %d %d %d %d %d %d %d \npanel tri %d %d %d %d %d %d %d %d %d \nend_surface\n\n';
fprintf(fid, words6, surfaces);
words7 = 'start_compartment cell%d \nsurface surf%d \npoint %d %d %d \nend_compartment \n\n';
fprintf(fid, words7, cells);
words8 = ' mol %d one u u u\n mol %d two u u u\n mol %d three u u u\n mol %d four u u u\n\n';
fprintf(fid, words8, no_mols, no_mols, no_mols, no_mols);
words9 = '\noutput_root results/ \noutput_files trajs1.txt trajs2.txt trajs3.txt trajs4.txt \ncmd n 1 listmols3 one trajs1.txt\ncmd n 1 listmols3 two trajs2.txt\ncmd n 1 listmols3 three trajs3.txt\ncmd n 1 listmols3 four trajs4.txt\n\nend_file\n';
fprintf(fid, words9);

% close the file
fclose (fid)

end

