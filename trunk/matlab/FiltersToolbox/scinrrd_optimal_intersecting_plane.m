function [ mopt, vopt, avals, mvals, vvals ] = ...
    scinrrd_optimal_intersecting_plane( nrrd, m0, v0, params )
% SCINRRD_OPTIMAL_INTERSECTING_PLANE  Optimise intersection plane for SCI
% NRRD segmentation mask
%
% [MOPT, VOPT, AVALS, MVALS, VVALS] = ...
%      SCINRRD_OPTIMAL_INTERSECTING_PLANE(NRRD, M0, V0, PARAMS)
%
%   This function computes the plane that intersects a SCI NRRD
%   segmentation mask in a way that minimizes the segmentation area
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
%   PARAMS is a struct with optimisation parameters:
%
%   * PARAMS.TYPE is a string. 'local' (default) means that the
%   intersected area will be minimised using multidimensional unconstrained
%   nonlinear minimization (Nelder-Mead, fminsearch). 'global' means that
%   area values will be systematically computed over a range of plane
%   inclination values.
%
%   * PARAMS.RAD can be a scalar or 2-vector:
%
%     scalar:   RAD is the radius of a 2D smoothing disk. The 2D image
%               obtained from intersecting the 3D volume by the plane at
%               each optimization iteration will be smoothed out using the
%               disk before computing the area.
%     2-vector: The smoothing element is a 3D ball. Instead of smoothing
%               a 2D image at each iteration, the whole 3D image volume is
%               smoothed with the ball at the beginning. RAD(1) is the ball
%               radius in the XY-plane. RAD(2) is the ball height.
%
%   * PARAMS.RANGE (global optimisation only): Angular range. The azimuth
%     of V0 will be changed RANGE(1) rad in any direction to compute
%     area values. The elevation of V0 will be changed RANGE(2) rad.
%     (Default RANGE(1) = RANGE(2) = 0.5236 rad = 30º).
%     
%   * PARAMS.N (global optimisation only): Number of azimuth (N(1)) or
%   elevation (N(2)) samples. (Default N(1) = N(2) = 61).
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
%   MVALS is a volume with the coordinates of M at the different tested
%   values. If optimisation is 'local', then MVALS is a matrix where each
%   column is the centroid at each iteration step. If optimisation is
%   'global', then MVALS is a volume where MVALS(T,P,:) is the centroid for
%   the normal vector VVALS(T,P,:).
%
%   VVALS is a volume like MVALS, only for the normal vectors.
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

%% Checks and initialization

% check arguments
error( nargchk( 2, 4, nargin, 'struct' ) );
error( nargoutchk( 0, 5, nargout, 'struct' ) );

% defaults
if ( nargin < 3 || isempty( v0 ) )
    v0 = [ 0 0 1 ];
end
if ( nargin < 4 || isempty( params ) )
    params.type = 'local';
    params.rad = [];
    params.range = [ 30 30 ] / 180 * pi;
    params.n = [ 61 61 ];
    se = [];
end
if (~isfield(params, 'type' ))
    params.type = 'local';
end
if (~isfield(params, 'rad' ))
    params.rad = [];
    se = [];
end
if (~isfield(params, 'range' ))
    params.range = [ 30 30 ] / 180 * pi;
end
if (~isfield(params, 'n' ))
    params.n = [ 61 61 ];
end

% prevent user entering rotation matrix instead of initial vector by
% mistake
if ( size( v0, 2 ) ~= 1 || size( v0, 1 ) ~= 3 )
    error( 'V0 must be a column 3-vector' )
end

% remove the dummy dimension and convert image data to double
nrrd = scinrrd_squeeze( nrrd, true );

% convert radius from real world size into number of pixels
params.rad = round( ...
    params.rad ./ [ nrrd.axis( 1:length( params.rad ) ).spacing ] );

% create disk for dilation/erosion if we are going to smooth in 2D
if ( length( params.rad ) == 1 )
    se = strel( 'disk', params.rad );
end

% 3D smoothing of the segmentation edges
if ( length( params.rad ) == 2 )
    se = strel( 'ball', params.rad(1), params.rad(2) );
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

% convert Cartesian coordinates into spherical coordinates (length has to
% be one); note: we use phi for azimuth, and theta for elevation, contrary
% to Matlab's naming convention
[ phi0, theta0 ] = cart2sph( v0(1), v0(2), v0(3) );

% group spherical coordinates into vector
alpha0 = [ phi0, theta0 ];

% init variables to keep track of the evolution of area values in the
% optimisation
avals = [];
mvals = [];
vvals = [];
    
% local or global optimisation
if strcmp( params.type, 'local' )

    % run optimisation to find minimum area; note that v0 is the only
    % variable optimised, but the rest (nrrd, x, y, z, m, rad, se) are
    % available to segmented_area_of_intersection() because the latter is a
    % subfunction
    alpha = fminsearch(@segmented_area_of_intersection, alpha0);
    
    % convert result from spherical to Carterian coordinates
    [ aux1 aux2 aux3 ] = sph2cart( alpha(1), alpha(2), 1.0 );
    vopt = [ aux1 aux2 aux3 ];
    
    % final centroid of the intersecting plane
    mopt = mvals( :, end );

elseif strcmp( params.type, 'global' )
    
    % interval of azimuth angle values, phi \in [-180º, 180º] or \in 
    % [0, 360º]
    phimin = phi0 - abs( params.range(1) );
    phimax = phi0 + abs( params.range(1) );
    
    % interval of elevation angle values, theta \in [-90º, 90º]
    thmin = max( theta0 - abs( params.range(2) ), -pi/2 );
    thmax = min( theta0 + abs( params.range(2) ), pi/2 );
    
    % sample angle intervals
    phivals = linspace( phimin, phimax, params.n(1) );
    thetavals = linspace( thmin, thmax, params.n(2) );
    
    % create matrices to save outputs; note that for each area value we
    % need to save a 3-vector with the rotation point, and a 3-vector with
    % the normal plane
    avals = zeros( length(thetavals), length(phivals) );
    mvals = zeros( length(thetavals), length(phivals), 3 );
    vvals = zeros( length(thetavals), length(phivals), 3 );
    
    % compute area for each combination of elevation and azimuth angles
    for T = 1:length( thetavals ) % elevation
        for P = 1:length( phivals ) % azimuth
            [ a, mnew, v ] = segmented_area_of_intersection( ...
                [ phivals(P) thetavals(T) ] );
            
            % put values in output matrices
            avals(T, P) = a;
            mvals(T, P, :) = mnew;
            vvals(T, P, :) = v;
            
        end
    end
    
    % find minimum area
    [foo, idx] = min( avals(:) );
    
    % convert linear index to multiple subscripts
    [T, P] = ind2sub( size(avals), idx );
    
    % output optimal plane
    mopt = squeeze( mvals(T, P, :) );
    vopt = squeeze( vvals(T, P, :) );
    
else
    error( [ 'Optimisation type not implemented: ' params.type ] )
end

    %% Objective function (the function we are trying to minimise)
    
    % rotate plane, intersect with image, and compute segmented area
    function [a, mnew, v] = segmented_area_of_intersection(alpha)
        
        % convert spherical to Carterian coordinates
        [ aux1 aux2 aux3 ] = sph2cart( alpha(1), alpha(2), 1.0 );
        v = [ aux1 aux2 aux3 ]';
        
        % vector cannot be zero
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
        if ( length( params.rad ) == 1 )
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
        rotmat = vec2rotmat( v );
  
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
        mnew = (mnew' + m0)';
        
        % for the global algorithm, values are recorded in a different way
        if strcmp( params.type, 'local' )
            % kept track of optimisation evolution
            avals = [ avals a ];
            vvals = [ vvals v ];
            mvals = [ mvals mnew ];
        end
        
    end
end
