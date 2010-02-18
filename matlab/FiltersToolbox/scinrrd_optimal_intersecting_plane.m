function [ mopt, vopt, avals, mvals, vvals ] = scinrrd_optimal_intersecting_plane( nrrd, m0, v0, rad )
% SCINRRD_OPTIMAL_INTERSECTING_PLANE  Optimise intersection plane for SCI
% NRRD segmentation mask
%
% [MOPT, VOPT, AVALS, MVALS, VVALS] = SCINRRD_OPTIMAL_INTERSECTING_PLANE(NRRD, M0, V0, RAD)
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
%   M0 is the rotation centroid. This centroid will not change.
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
%   MOPT is the centroid of the optimal intersection.
%
%   VOPT is the normal vector to the optimal plane.
%
%   Note that because M0 and MOPT are both contained in the optimal plane,
%   the duples (M0, VOPT) and (MOPT, VOPT) define the same optimal plane.
%
%   AVALS is a vector with the record of area values from the optimization.
%
%   MVALS is a matrix with the coordinates of MOPT at each optimization
%   iteration.
%
%   VVALS is a matrix with the record of VOPT at each optimization
%   iteration.
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

% Author(s): Ramon Casero, Vicente Grau
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

%% Checks and initialization

% check arguments
error( nargchk( 2, 4, nargin, 'struct' ) );
error( nargoutchk( 0, 5, nargout, 'struct' ) );

% defaults
if ( nargin < 3 || isempty( v0 ) )
    v0 = [ 0 0 1 ];
end
if ( nargin < 4 || isempty( rad ) )
    rad = [];
    se = [];
end

% prevent user entering rotation matrix instead of initial vector
if ( size( v0, 2 ) ~= 1 || size( v0, 1 ) ~= 3 )
    error( 'V0 must be a column 3-vector' )
end

% remove the dummy dimension and convert image data to double
nrrd = scinrrd_squeeze( nrrd, true );

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

% generate 3D grid of coordinates
[ x, y, z ] = scinrrd_ndgrid( nrrd );

% % DEBUG: compute intersection of NRRD volume with the initial plane
% % (if you want to visualize the image as in Seg3D, you need to do 'axis
% % xy')
% im = scinrrd_intersect_plane(nrrd, m0, v0, x, y, z);

%% Optimisation of the intersection area

% init variables to keep track of the evolution of area values in the
% optimisation
avals = [];
mvals = [];
vvals = [];

% convert Cartesian coordinates into spherical coordinates (length has to
% be one); note: we use phi for azimuth, and theta for elevation, contrary
% to Matlab's convention
[ phi0, theta0 ] = cart2sph( v0(1), v0(2), v0(3) );

% group spherical coordinates into vector
alpha0 = [ phi0, theta0 ];

% run optimisation to find minimum area; note that v0 is the only variable
% optimised, but the rest (nrrd, x, y, z, m, rad, se) are available to
% segmented_area_of_interest() because the latter is a subfunction
alpha = fminsearch(@segmented_area_of_intersection, alpha0);

% convert result from spherical to Carterian coordinates
[ aux1 aux2 aux3 ] = sph2cart( alpha(1), alpha(2), 1.0 );
vopt = [ aux1 aux2 aux3 ];

% final centroid of the intersecting plane
mopt = mvals( :, end );

    %% Objective function (the function we are trying to minimise)
    
    % rotate plane, intersect with image, and compute segmented area
    function a = segmented_area_of_intersection(alpha)
        
        % convert spherical to Carterian coordinates
        [ aux1 aux2 aux3 ] = sph2cart( alpha(1), alpha(2), 1.0 );
        v = [ aux1 aux2 aux3 ];
        
        % normalize vector
        if ( norm(v) == 0 )
            error( 'Normal vector to plane cannot be (0,0,0)' )
        end
        
        % this function cannot deal with vertical planes, because of a
        % singularity
        if ( v(3) == 0 )
            error( 'Intersecting plane cannot be vertical' )
        end

        % compute intersection of plane with volume
        [ im, zp, xp, yp ] = scinrrd_intersect_plane(nrrd, m0, v, x, y, z);
        
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
        
%         % DEBUG: visualize intersection projected onto horizontal plane
%         hold off
%         imagesc(xp(:), yp(:), im > 0)
%         hold on
%         % DEBUG: compute and plot convex hull
%         idx2 = convhull( xps, yps );
%         vxs = xps(idx2);
%         vys = yps(idx2);
%         plot(vxs, vys, 'w')
%         xlabel( 'x (m)' )
%         ylabel( 'y (m)' )
%         pause
        
        % compute a rotation matrix from the Cartesian system to the
        % rotated plane
        rotmat = vec2rotmat( v' );
  
        % we are now seeing the rotated plane projected onto the horizontal
        % plane, i.e. we see the segmentation mask in perspective.
        % In order to see the true area of the segmentation mask, we need
        % to change the system of coordinates so that the rotated plane
        % becames the XY plane
        
        % first, move segmented points so that centroid is at (0,0,0)...
        xps = xps - m0(1);
        yps = yps - m0(2);
        zps = zps - m0(3);
        
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
%         plot(xps + m0(1), yps + m0(2), 'r*')
%         axis ij
        
        % compute convex hull (reuse idx2): note convex hull coordinates
        % are on projected space
        idx2 = convhull( xps, yps );
        vxs = xps(idx2);
        vys = yps(idx2);
        
        % compute x-,y-coordinates centroid and area of polygon
        [ mnew, a ] = polycenter( vxs, vys );
        mnew(3) = 0;
        
        % the centroid is now on projected coordinates, but we need to put
        % it back on the real world coordinates
        mnew = rotmat * mnew';
        mnew = mnew' + m0;
        
        % kept track of optimisation evolution
        avals = [ avals a ];
        vvals = [ vvals v' ];
        mvals = [ mvals mnew' ];
        
    end
end
