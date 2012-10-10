function [ vopt, mopt, avals, vvals, mvals ] = scinrrd_optimal_intersecting_plane( nrrd, z0, v0, rad )
% Note: Old implementation of function
% scinrrd_optimal_intersecting_plane(). Backup while we transition to
% a definite method.
%
% SCINRRD_OPTIMAL_INTERSECTING_PLANE  Optimise intersection plane for SCI
% NRRD segmentation mask
%
% [VOPT, MOPT, AVALS, VVALS, MVALS] = SCINRRD_OPTIMAL_INTERSECTING_PLANE(NRRD, Z0, V0, RAD)
%
%   This function computes the plane that intersects a SCI NRRD
%   segmentation mask in a way that minimizers the segmentation area
%   intersected by the plane. That is, in some sense in finds the plane
%   more orthogonal to the segmented volume.
%
%   (Note that the area is computed on the convex hull of the plane
%   intersection with the volume.)
%
%   NRRD is the SCI NRRD struct.
%
%   Z0 is a z-coordinate value. The rotation centroid for the plane will be
%   at Z0 height.
%
%   V0 is a 3-vector that describes the normal vector to the initial
%   intersecting plane. By default, the initial plane is horizontal.
%
%   RAD can be a scalar or 2-vector:
%     scalar:   RAD is the radius of a 2D smoothing disk. The 2D image
%               obtained from intersecting the 3D volume by the plane at
%               each optimization iteration will be smoothed out using the
%               disk before computing the area.
%     2-vector: The smoothing element is a 3D ball. Instead of smoothing
%               a 2D image at each iteration, the whole 3D image volume is
%               smoothed with the ball at the beginning. RAD(1) is the ball
%               radius in the XY-plane. RAD(2) is the ball height.
%
%   VOPT is the normal vector that describes the plane at centroid MOPT so
%   that the intersection area is minimised.
%
%   AVALS is a vector with the record of area values from the optimization.
%
%   VVALS is a matrix with the record of normal vectors to the plane from
%   the optimization.
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

% Author(s): Ramon Casero <rcasero@gmail.com>, Vicente Grau
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
error( nargchk( 2, 4, nargin, 'struct' ) );
error( nargoutchk( 0, 5, nargout, 'struct' ) );

% defaults
if ( nargin < 3 || isempty( v0 ) )
    v0 = [ 0 0 1 ];
end
if ( nargin < 4 || isempty( rad ) )
    rad = [];
end

% remove the dummy dimension and convert image data to double
nrrd = scinrrd_squeeze( nrrd, true );

% compute image size
sz = size( nrrd.data );

% convert radius from real world size into number of pixels
rad = round( rad ./ [ nrrd.axis( 1:length( rad ) ).spacing ] );

% create disk for dilation/erosion if we are going to smooth in 2D
if ( length( rad ) == 1 )
    se = strel( 'disk', rad );
end

% 3D smoothing of the segmentation edges
if ( length( rad ) == 2 )
    se = strel( 'ball', rad(1), rad(2) );
    nrrd.data = imdilate( nrrd.data, se );
    nrrd.data = imerode( nrrd.data, se );
end

% get linear indices of segmented voxels
idx = find( nrrd.data );

% convert the linear indices to volume indices
[ix, iy, iz] = ind2sub( sz( 1:end ), idx );

% compute real world coordinates for those indices
coords = scinrrd_index2world( [ ix, iy, iz ], nrrd.axis );

% get tight frame around segmentation
cmin = min( coords );
cmax = max( coords );

% compute x-, y-coordinates of centroid
m = mean( coords );

% defaults
% 2/3 from the bottom (we use this as reference plane for the left
% ventricle in the heart)
if ( nargin < 2 || isempty( z0 ) )
    z0 = cmax(3) - ( cmax(3) - cmin(3) ) / 3;
end
m(3) = z0;

% generate 3D grid of coordinates
[ x, y, z ] = scinrrd_ndgrid( nrrd );

% compute intersection of NRRD volume with the hrozontal 2D plane at height
% z0
% (if you want to visualize the image as in Seg3D, you need to do 'axis
% xy')
[ im, zp0, xp0, yp0 ] = scinrrd_intersect_plane(nrrd, m, [0 0 1], x, y, z);

% get linear indices of segmented voxels in the 2D intersection
idx2 = find( im );

% get coordinates of segmented voxels
xp = xp0( idx2 );
yp = yp0( idx2 );
zp = zp0( idx2 );

% % DEBUG: to visualize segmentation mask in real world coordinates
% hold off
% imagesc(xp0(:), yp0(:), im)
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

% Note: zp0 is only used for debugging
[ vopt, mopt, avals, vvals, mvals ] = optimise_plane_rotation(v0, nrrd, x, y, z, m, rad, se);

end % function scirunnrrd_intersect_plane

% function to optimise the rotation of the plane so that it minimises the
% segmented area
function [ v, m, avals, vvals, mvals ] = optimise_plane_rotation(v0, nrrd, x, y, z, m, rad, se)

% init variables to keep track of the evolution of area values in the
% optimisation
avals = [];
vvals = [];
mvals = [];

% we need to do this for the first centroid update (the value of m will be
% updated from within the optimization function)
mnew = m;

% run optimisation to find minimum area
v = fminsearch(@segmented_area_of_intersection, v0);

% normalize vector
v = v / norm( v );

% rotate plane, intersect with image, and compute segmented area
    function a = segmented_area_of_intersection(v)

        % this function cannot deal with vertical planes, because of a
        % singularity
        if ( v(3) == 0 )
            error( 'Intersecting plane cannot be vertical' )
        end

        % update the centroid
        m = mnew;
        
        % normalize vector
        if ( norm(v) == 0 )
            error( 'Normal vector to plane cannot be (0,0,0)' )
        end
        v = v / norm(v);
        
        % compute intersection of plane with volume
        [ im, zp ] = scinrrd_intersect_plane(nrrd, m, v, x, y, z);
        xp = x( :, :, 1 );
        yp = y( :, :, 1 );
        
        % 2D smoothing of the segmentation edges
        if ( length( rad ) == 1 )
            im = imdilate( im, se );
            im = imerode( im, se );
        end
        
%         % DEBUG: plot rotated plane
%         hold off
%         plot3( xp(:), yp(:), zp(:), '.r' )
        
        % find segmented voxels in the 2D cut
        idx = find( im );

        % get coordinates of segmented voxels
        xps = xp( idx );
        yps = yp( idx );
        zps = zp( idx );
        
        % DEBUG: to visualize segmentation mask in real world coordinates
        hold off
        imagesc(xp(:), yp(:), im > 0)
        hold on
%         plot(xps, yps, 'w*')
        plot(m(1), m(2), 'wo')
        % DEBUG: compute convex hull
        idx2 = convhull( xps, yps );
        vxs = xps(idx2);
        vys = yps(idx2);
        % DEBUG: plot convex hull
        plot(vxs, vys, 'w')
        xlabel( 'x (m)' )
        ylabel( 'y (m)' )
        pause
        
        
        % compute a rotation matrix that will transform the Cartesian
        % system of reference onto another one where the XY plane is the
        % rotated plane; note that to compute the area, we can use any pair
        % of orthogonal vectors on the rotated plane, plus the unique
        % vector orthogonal to the plane
        
        % compute one arbitrary vector contained by the plane, e.g. the
        % corresponding to coordinates [1, 0, ...]. By definition, any
        % vector contained in the plane has to be orthogonal to v
        v2 = [ 1, 0, ...
            v(1)/v(3)*(m(1)-1) + v(2)/v(3)*(m(2)-0) + m(3) ];

        % normalize the vector
        v2 = v2 / norm( v2 );
        
        % make the third vector orthogonal to the previous two
        v3 = cross( v, v2 );
        
        % normalize vector to correct numerical errors
        v3 = v3 / norm( v3 );
        
        % the rotation matrix is the rotated orthonormal basis, so that if
        % we do rotmat * eye(3), the Cartesian basis becomes the rotated
        % basis
        rotmat = [ v2'  v3'  v' ];
        
        % we are now seeing the rotated plane projected onto the horizontal
        % plane, i.e. we see the segmentation mask in perspective.
        % In order to see the true area of the segmentation mask, we need
        % to change the system of coordinates so that the rotated plane
        % becames the XY plane
        
        % first, move segmented points so that centroid is at (0,0,0)...
        xps = xps - m(1);
        yps = yps - m(2);
        zps = zps - m(3);
        
        % ...second, make the rotated plane horizontal, by inverting the
        % rotation...
        xyzps = [ xps(:) yps(:) zps(:) ] * rotmat;
        xps = xyzps( :, 1 );
        yps = xyzps( :, 2 );
        zps = xyzps( :, 3 );
        
        % if everything has gone alright, then the z-coordinate of xyzp
        % should be zero (+numerical errors), because the rotated plane is
        % now the XY plane
        assert( abs( min( zps ) ) < 1e-10 )
        assert( abs( max( zps ) ) < 1e-10 )
        
%         % DEBUG: visualize segmentation mask in real world coordinates
%         hold off
%         imagesc(xp(:), yp(:), im > 0)
%         hold on
%         plot(xps + m(1), yps + m(2), 'w*')
        
        % compute convex hull (reuse idx2): note convex hull coordinates
        % are on projected space, image is in original space
        idx2 = convhull( xps, yps );
        vxs = xps(idx2);
        vys = yps(idx2);
        
        % compute x-,y-coordinates centroid and area of polygon
        [ mnew, a ] = polycenter( vxs, vys );
        mnew(3) = 0;
        
        % the centroid is now on projected coordinates, but we need to put
        % it back on the real world coordinates
        mnew = rotmat * mnew';
        mnew = mnew' + m;
        
        % the new centroid is now on the intersecting plane, generally at a
        % height ~= z0. Now we project the new centroid on the horizontal
        % plane at height == z0
        mnew(3) = m(3);
        
        % kept track of optimisation evolution
        avals = [ avals a ];
        vvals = [ vvals v' ];
        mvals = [ mvals mnew' ];
        
    end
end
