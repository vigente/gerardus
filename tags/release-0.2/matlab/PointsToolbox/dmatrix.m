function d = dmatrix( x, y )
% DMATRIX  Matrix of distances between vectors
%
%    D = DMATRIX(X)
%
%      D(i,j) is the Euclidean distance between X(:,i) and X(:,j)
%
%    D = DMATRIX(X,Y)
%
%      D(i,j) is the Euclidean distance between X(:,i) and Y(:,j)

% Author: Ramon Casero <rcasero@robots.ox.ac.uk>
% Copyright © 2007 Ramon Casero, Wolfson Medical Vision Laboratory, 
%                  University of Oxford
% v2.0 - 6 Sep 2007

% check arguments
msg = nargchk( 1, 2, nargin );
if ( ~isempty( msg ) ), error( msg ), end
msg = nargoutchk( 0, 1, nargout );
if ( ~isempty( msg ) ), error( msg ), end

%    D = DMATRIX(X)
if ( nargin == 1 )
    
    M = size( x, 2 );
    d = zeros( M );
    for I = 1:M
        for J = (I + 1):M
            d( I, J ) = sqrt( sum( ( x( :, I ) - x( :, J ) ) .^ 2 ) );
            d( J, I ) = d( I, J );
        end
    end

%    D = DMATRIX(X,Y)
else
    
    M = size( x, 2 );
    N = size( y, 2 );
    d = zeros( M, N );
    for I = 1:M
        for J = 1:N
            d( I, J ) = sqrt( sum( ( x( :, I ) - y( :, J ) ) .^ 2 ) );
        end
    end

end
