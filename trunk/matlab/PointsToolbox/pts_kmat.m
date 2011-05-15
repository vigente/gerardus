function k = pts_kmat( kernel, x, y )
% PTS_KMAT  Compute kernel matrix for Linear and Kernel Principal Component
% Analysis (PCA and KPCA)
%
% K = PTS_KMAT(KERNEL, X)
%
%    Compute kernel matrix for vectors X (a matrix of column vectors) using
%    KERNEL.
%
% K = PTS_KMAT(KERNEL, X, Y)
%
%    Compute kernel matrix for vectors X, Y (matrices of column vectors)
%    using KERNEL.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2009 University of Oxford
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
error( nargchk( 2, 3, nargin, 'struct' ) );
error( nargoutchk( 0, 1, nargout, 'struct' ) );

switch kernel.type
    case 'polyh'
        q = kernel.params; % polynomial degree
        if ( nargin < 3 )
            k = ( x' * x ) .^ q;
        else
            k = ( x' * y ) .^ q;
        end
    case 'gauss'
        sigma = kernel.params; % standard deviation
        beta = - 1 / ( 2 * sigma^2 );
        if ( nargin < 3 )

            M = size( x, 2 ); % number of input vectors
            % init kernel matrix
            k = zeros( M );
            for I = 1:M
                for J = I:M
                    % compute upper half of kernel matrix
                    k( I, J ) = exp( beta * ...
                        sum( ( x( :, I ) - x( :, J ) ) .^ 2 ) );
                    % use symmetry for lower half
                    k( J, I ) = k( I, J );
                end
            end
            
        else
            
            M1 = size( x, 2 ); % number of input vectors
            M2 = size( y, 2 ); % number of input vectors
            % init kernel matrix
            k = zeros( M1, M2 );
            for I = 1:M1
                for J = 1:M2
                    % compute kernel matrix
                    k( I, J ) = exp( beta * ...
                        sum( ( x( :, I ) - y( :, J ) ) .^ 2 ) );
                end
            end

        end
        
        % check that we are not losing information in k due to lack of
        % numerical precision: the smallest number that Matlab can
        % represent is eps, and this corresponds to a squared distance of
        % -log(eps) = 36.0437 (in my machine); for larger distances, k=0
        % and then d^2 = Inf, i.e. the true value of d^2 is lost
        %
        % the solution is to make sigma bigger, so that the term in the
        % exp() is smaller in absolute value
        if any( k( : ) == 0 )
            warning( 'Values of K equal to 0, you need to make SIGMA bigger' )
        end

    otherwise
        error( 'Invalid type of kernel' )
end
