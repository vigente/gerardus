function [ qopt, mopt ] = scirunnrrd_intersect_plane( im, z0 )
% SCIRUNNRRD_INTERSECT_PLANE  Optimise intersection plane for SCI NRRD
% segmentation mask
%
% [QOPT, MOPT] = SCIRUNNRRD_INTERSECT_PLANE(IM, Z0)
%
%   This function computes the plane that intersects a SCI NRRD
%   segmentation mask in a way that minimizers the segmentation area
%   intersected by the plane. That is, in some sense in finds the plane
%   more orthogonal to the segmented volume.
%
%   (Note that the area is computed on the convex hull of the plane
%   intersection with the volume.)
%
%   IM is the SCI NRRD struct.
%
%   Z0 is a z-coordinate value. The rotation centroid for the plane will be
%   at Z0 height.
%
%   QOPT is the quaternion that describes the optimal rotation of the plane
%   around centroid MOPT so that the intersection area is minimised.
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

% Copyright © 2010 University of Oxford
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

% get image vertices
vmin = [ im.axis.min ];
vmax = [ im.axis.max ];

% get voxel size
res = [ im.axis.spacing ];

% squeeze the non-used first dimension of data
im.data = double( squeeze( im.data ) );
vmin = vmin( 2:end );
vmax = vmax( 2:end );
res = res( 2:end );

% compute image size
sz = size( im.data );

% get linear indices of segmented voxels
idx = find( im.data );

% convert the linear indices to volume indices
[ix, iy, iz] = ind2sub( sz( 2:end ), idx );

% compute real world coordinates for those indices
coords = scirunnrrd_index2world( [ ix, iy, iz ], im.axis );

% get tight frame around segmentation
cmin = min( coords );
cmax = max( coords );

% compute centroid
m = mean( coords );

% defaults
% 2/3 from the bottom (we use this as reference plane for the left
% ventricle in the heart)
if ( nargin < 2 || isempty( z0 ) )
    z0 = cmax(3) - ( cmax(3) - cmin(3) ) / 3;
end

% generate 3D grid of coordinates
% note the inversion of x and y so that x values change with columns, and y
% values change with rows
[ y, x, z ] = ndgrid( vmin(1):res(1):vmax(1), ...
    vmin(2):res(2):vmax(2), vmin(3):res(3):vmax(3) );

% compute a horizontal 2D plane at the height z0
[ yp0, xp0, zp0 ] = ndgrid( vmin(1):res(1):vmax(1), ...
    vmin(2):res(2):vmax(2), z0 );

% compute intersection of image with the 2D plane
% (if you want to visualize the image as in Seg3D, you need to do 'axis
% xy')
% note: this is counterintuitive, because it seems that the first
% coordinate (rows) should be sampled by y. But the way x and y are
% defined, this is how it works
im2 = interpn( x, y, z, im.data, xp0, yp0, zp0, ...
    'linear' );

% get linear indices of segmented voxels in the 2D intersection
idx2 = find( im2 );

% get coordinates of segmented voxels
xp = xp0( idx2 );
yp = yp0( idx2 );
zp = zp0( idx2 );

% % DEBUG: to visualize segmentation mask in real world coordinates
% hold off
% imagesc(xp0(:), yp0(:), im2)
% hold on
% plot(xp, yp, 'w*')

% compute convex hull (reuse idx2)
idx2 = convhull( xp, yp );
vx = xp(idx2);
vy = yp(idx2);

% % DEBUG: plot convex hull
% plot(vx, vy, 'r')

% compute centroid and area of polygon
[ m, a ] = polycenter( vx, vy );

% % DEBUG: plot centroid
% plot(m(1), m(2), 'ko')

% place centroid at correct height
m(3) = z0;

% initialize quaternion with no rotation
q0 = [ 1 0 0 0 ];

% Note: zp0 is only used for debugging
[ qopt, mopt ] = optimise_plane_rotation(q0, im, x, y, z, m, xp0, yp0, zp0);

end % function scirunnrrd_intersect_plane

% function to optimise the rotation of the plane so that it minimises the
% segmented area
function [q, m] = optimise_plane_rotation(q0, im, x, y, z, m, xp0, yp0, zp0)

% we need to do this for the first centroid update
mnew = m;

% run optimisation to find minimum area
[q,a] = lsqnonlin(@segmented_area, q0);

% rotate plane, intersect with image, and compute segmented area
    function a = segmented_area(q)
        
        % update the centroid
        m = mnew;
        
        % compute image size
        sz = size( im.data );
        
        % normalize quaternion
        if ( norm(q) ~= 0 )
            q = q / norm(q);
        end
        
        % compute rotation matrix for quaternion
        rotmat = quaternion2matrix( q );
        
        % A plane with orthogonal vector n that goes through the centroid m
        % can be expressed by the formula
        %
        %   nx(x-mx) + ny(y-my) + nz(z-mz) = 0
        %
        % We know the horizontal limits of the image. So in order to
        % compute the intersection height of the rotated plane at each one
        % of 2 vertices that delimitate the image horizontally, we use the
        % derived expression
        %
        %    z = nx/nz(mx-x) + ny/nz(my-y) + mz
        %
        % This expression is invalid when nz=0, i.e. the plane has been
        % rotated to make it "vertical"
        
        % the normal vector associated to a horizontal plane is vertical
        n = [ 0 0 1 ]';
        
        % rotate the normal vector (note that by definition it's at the
        % centroid already)
        n = rotmat * n;
        
        if ( n(3) == 0 )
            error( 'Rotation makes plane vertical' )
        end
        
        % the horizontal grid for the rotated plane will be the same as for
        % the horizontal plane
        xp = xp0;
        yp = yp0;
        
        % compute the corresponding heights at the horizontal grid, i.e.
        % the rotated plane
        zp = n(1)/n(3)*(m(1)-xp0) + n(2)/n(3)*(m(2)-yp0) + m(3);
        
        %         % DEBUG: plot horizontal plane
        %         hold off
        %         plot3( xp0(:), yp0(:), zp0(:), '.' )
        %         % DEBUG: plot rotated plane
        %         hold on
        %         plot3( xp(:), yp(:), zp(:), '.r' )
        
        % reshape coordinates as matrix, otherwise sampling is not in the
        % right order
        xp = reshape( xp, sz(1:2) );
        yp = reshape( yp, sz(1:2) );
        zp = reshape( zp, sz(1:2) );
        
        % compute intersection of image with the 2D plane
        % (if you want to visualize the image as in Seg3D, you need to do
        % 'axis xy')
        % note: this is counterintuitive, because it seems that the first
        % coordinate (rows) should be sampled by y. But the way x and y are
        % defined, this is how it works
        im2 = interpn( x, y, z, im.data, xp, yp, zp, ...
            'linear', 0 );
        
        % find segmented voxels in the 2D cut
        idx2 = find( im2 );
        
        % get coordinates of segmented voxels
        xp = xp0( idx2 );
        yp = yp0( idx2 );
        zp = zp( idx2 );
        
        %         % DEBUG: to visualize segmentation mask in real world coordinates
        %         hold off
        %         imagesc(xp0(:), yp0(:), im2)
        %         hold on
        %         plot(xp, yp, 'w*')
        %         plot(m(1), m(2), 'ko')
        
        % we are now seeing the rotated plane projected onto the horizontal
        % plane, i.e. we see the segmentation mask in perspective.
        % In order to see the true area of the segmentation mask, we need
        % to change the system of coordinates so that (0,0,0) is on the
        % rotation centroid, and the rotated plane becames the XY plane
        
        % first, move plane so that centroid is at (0,0,0)...
        xp = xp - m(1);
        yp = yp - m(2);
        zp = zp - m(3);
        
        % ...second, make the rotated plane horizontal, by inverting the
        % rotation...
        xyzp = [ xp(:) yp(:) zp(:) ] * rotmat;
        xp = xyzp( :, 1 );
        yp = xyzp( :, 2 );
        zp = xyzp( :, 3 );
        
        % if everything has gone alright, then the z-coordinate of xyzp
        % should be zero (+numerical errors), because the rotated plane is
        % now the XY plane
        assert( abs( min( zp ) ) < 1e-10 )
        assert( abs( max( zp ) ) < 1e-10 )
        
        % ... and finally, move the plane back to the original centroid
        xp = xp + m(1);
        yp = yp + m(2);
        zp = m(3);
        
        %         % DEBUG: visualize segmentation mask in real world coordinates
        %         hold off
        %         imagesc(xp0(:), yp0(:), im2)
        %         hold on
        %         plot(xp, yp, 'w*')
        %         pause
        
        % compute convex hull (reuse idx2)
        idx2 = convhull( xp, yp );
        vx = xp(idx2);
        vy = yp(idx2);
        
        % % DEBUG: plot convex hull
        % plot(vx, vy, 'r')
        
        % compute x-,y-coordinates centroid and area of polygon
        [ mnew, a ] = polycenter( vx, vy );
        mnew(3) = m(3);
        
        % the centroid is now on projected coordinates, but we need to put
        % it back on the real world coordinates
        mnew = mnew - m;
        mnew = rotmat * mnew';
        mnew = mnew' + m;
        
        % we also need to ignore the z-coordinate of the new centroid
        mnew(3) = m(3);
        
    end
end

