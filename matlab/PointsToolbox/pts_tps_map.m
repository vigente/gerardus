function [y, w] = pts_tps_map( s, t, x, w, FAST, PROGRESS )
% PTS_TPS_MAP  Interpolate/warp/map N-dimensional points using a thin-plate
% spline transformation
%
% [Y, W] = PTS_TPS_MAP(S, T, X, W)
%
%    S is a (P,Ds,N)-volume where each (:,:,i)-matrix has the coordinates
%    of the source points that define the warp.
%
%    T is a (P,Dt,N)-volume where each (:,:,i)-matrix has the coordinates
%    of the target points that define the warp.
%
%      (P is the number of points, Ds and Dt are the dimension and N is the
%      number of configurations).
%
%    X is a (P,Ds,N)-volume where each (:,:,i)-matrix has the coordinates
%    of the points to be interpolated.
%
%    Y is a (P,Dt,N)-volume with the interpolated points.
%
%    X can also be a struct with 3 volumes: X.in (P1,Ds,N), X.endo (P2,Ds,N),
%    X.epi (P3,Ds,N). In this case, 3 warps are computed: With all points
%    for X.in; with the first half of S and T points for X.endo; and with
%    the second half of S and T points for X.epi.
%
%    W is a D2-column matrix where each column has the weight and affine
%    weights vector computed with PTS_TSP_WEIGHTS (if W is empty or
%    missing, then W is computed internally).
%
%    If X is a struct, then W has to be a similar struct with the weights
%    for each of the warps.
%
% ... = PTS_TPS_MAP(S, T, X, W, FAST)
%
%    To speed up computations, all points in the same configuration are
%    warped together using vector and matrix operations.
%
%    This is more memory expensive. If for some reason memory becomes a
%    problem, using FAST=false it is possible to run the algorithm using
%    a "for" loop to warp 1 point at a time. This is roughly 14 times
%    slower. By default the faster method is used, FAST=true.
%
% ... = PTS_TPS_MAP(S, T, X, W, FAST, PROGRESS)
%
%    Making PROGRESS=true, and if you have function STATUSBAR() by
%    Leutenegger Marcel, you will see a progress bar with an estimate of
%    the time remaining to completion if you are using the slow method.
%    
%
% See also: pts_tps_weights.

% Author: Ramon Casero <rcasero@gmail.com>
% v6.3 - 25 Jan 2011
% Copyright © 2006-2011 University of Oxford
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
error( nargchk( 3, 6, nargin ) );
error( nargoutchk( 0, 2, nargout ) );

% defaults
if ( nargin < 5 || isempty( FAST ) )
    FAST = true;
end
if ( nargin < 6 || isempty( PROGRESS ) )
    PROGRESS = false;
end

% struct and sizes check
P = size( s, 1 ); % number of control points
if ( P ~= size( t, 1 ) )
    error( 'S and T must have the same number of points' )
end

% compute weights if needed
N = size( s, 3 ); % number of frames
if ( nargin < 4 || isempty( w ) )
    if isstruct( x ) % triple warp
        
        for I = 1:N
            w( I ).in = pts_tps_weights( s( :, :, I ), t( :, :, I ) );
            w( I ).endo = pts_tps_weights( s( 1:(P/2), :, I ), ...
                t( 1:(P/2), :, I ) );
            w( I ).epi = pts_tps_weights( s( (P/2+1):end, :, I ), ...
                t( (P/2+1):end, :, I ) );
        end
        
    else
        w = pts_tps_weights( s, t );
    end
end

% deal with empty point sets
if isempty( x )
    y = [];
    return;
end

if isstruct( x ) % triple warp
    if any( ~isfield( x, {'endo', 'epi', 'in'} ) )
        error( 'Invalid fields for struct in X' )
    end
    
    % compute each warp separately
    y( N ) = struct( ...
        'in', zeros( size( x( 1 ).in ) ), ...
        'endo', zeros( size( x( 1 ).endo ) ), ...
        'epi', zeros( size( x( 1 ).epi ) ) );
    for I = 1:N
        y( I ).in = pts_tps_map( s( :, :, I ), t( :, :, I ), ...
            x( I ).in, w( I ).in, FAST );
        y( I ).endo = pts_tps_map( s( 1:(P/2), :, I ), ...
            t( 1:(P/2), :, I ), x( I ).endo, w( I ).endo, FAST );
        y( I ).epi = pts_tps_map( s( (P/2+1):end, :, I ), ...
            t( (P/2+1):end, :, I ), x( I ).epi, w( I ).epi, FAST );
    end
    
    return
end

% dimensionality of points
D = size( s, 2 );

if ( D ~= size( x, 2 ) )
    error( [ 'Dimensionality of S is ' num2str( D ) ', of X is ' num2str( size( x, 2 ) ) ] )
end
if ( N ~= size( t, 3 ) || N ~= size( x, 3 ) )
    error( 'S, T and X must have the same number of point configurations' )
end

% init output
y = zeros( size( x, 1 ), size( w, 2 ) );

PP = size( x, 1 ); % PP: number of points to be interpolated

if FAST

    % warp every input point
    for I = 1:N % loop frames

        % init aux matrix
        u = zeros( P * PP, D );

        % interleave points to be warped
        for J = 1:D
            u( :, J ) = reshape( repmat( x( :, J, I ), 1, P )', P * PP, 1 );
        end

        % compute norm(Pi - (x,y)).^2
        u = sum( ( u - repmat( s( :, :, I ), PP, 1 ) ) .^ 2, 2 );

        % reshape to have 1 row vector per point to be interpolated
        u = reshape( u, P, PP )';

        % thin-plate spline distance function
        % U(r) = r^2 log10(r)
        % note: it's faster to compute the log10 this way than directly
        warning( 'off', 'MATLAB:log:logOfZero' );
        u = 0.5 * u .* log( u ) * ( 1/log(10) );
        u( isnan( u ) ) = 0;
        warning( 'on', 'MATLAB:log:logOfZero' );

        % factor by weights f(x,y) = a1 + ax*x + ay*y + sum( wi*U(|Pi - (x,y)|) )
        y( :, :, I ) = [ u , ones( PP, 1 ) , x( :, :, I ) ] * w( :, :, I );

    end

else % loop throug each point; this is slower but less memory consuming
    
    % DEBUG
    error('Implementation of SLOW mode is broken and needs to be fixed. Use FAST mode')

    % show status bar
    if PROGRESS
        delete( statusbar )
        bar = statusbar( 'Progress bar...' );
    else
        bar = [];
    end
    
    % warp every input point
    for I = 1:N % loop frames
        for J = 1:PP
            % compute norm(Pi - (x,y)).^2
            % U(r) = r^2 log10(r)
            % note: it's faster to compute the log10 this way than directly
            u = sum( ...
                ( repmat( x( J, :, I ), P, 1 ) - s( :, :, I ) ) .^ 2, 2 )';
            warning( 'off', 'MATLAB:log:logOfZero' );
            u = 0.5 * u .* log( u ) * ( 1/log(10) );
            u( isnan( u ) ) = 0;
            warning( 'on', 'MATLAB:log:logOfZero' );

            % factor by weights f(x,y) = a1 + ax*x + ay*y + sum( wi*U(|Pi - (x,y)|) )
            y( J, :, I ) = ...
                [ u , ...
                1 , ...
                x( J, :, I ) ] ... % x, y, z, etc
                * w( :, :, I );

            % update progress bar every 10th of the total number of points
            if mod( J, round(PP/10) )
                if isempty( statusbar( ((I-1)*PP+J)/(N*PP), bar ) )
                    break
                end
            end

        end
        
%         % update progress bar every frame
%         if isempty( statusbar( I/N, bar ) )
%             break
%         end
        
    end
    
    % delete progress bar
    delete( bar );

end
