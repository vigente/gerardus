function [ v, d ] = pts_pca( x, kernel )
% PTS_PCA  Linear and Kernel Principal Component Analysis (PCA and KPCA)
%
% [ V, D ] = PTS_PCA(X)
%
%   Center X and compute eigenvectors V and eigenvalues D using linear PCA.
%
%   X is a matrix where each column is a vector. Each row is a sample of a
%   variable.
%
%   The eigenproblem solved is
%
%     S * v = d * v
%
%   where v is each eigenvector, and d each corresponding eigenvalue. S is
%   the biased covariance matrix
%
%     S = 1/M * Xc * Xc'
%
%   where Xc is the matrix of centered vectors X.
%
%   NOTE: If there are less vectors than variables, or kernel PCA is
%   selected, MDS is used to speed up computation of eigenvectors.
%
%   NOTE: Tiny eigenvalues can give numerical errors, so eigenvalues <
%   1e-13 and the corresponding eigenvectors are removed from the output.
%   Negative eigenvalues are assumed to arise from numerical errors, and
%   removed too.
%
%   NOTE: Matlab's function cov() computes the _unbiased_ covariance.
%
% [ A, DA ] = PTS_PCA(X, KERNEL)
%
%   Generalization to any type of kernel. KERNEL is a struct with the
%   kernel description
%
%     Linear PCA, as above
%
%       KERNEL.type = 'linear'
%
%     Kernel PCA, Polynomial homogeneous order q
%
%       KERNEL.type = 'polyh'
%       KERNEL.params = q
%
%     Kernel PCA, Gaussian with standard deviation sigma
%
%       KERNEL.type = 'gauss'
%       KERNEL.params = sigma
%
%   Because in general we want to avoid any explicit computations in
%   feature space, we solve the eigenproblem [1]
%
%     Kc * a = da * a
%
%   where a is each "coefficient eigenvector", and da = m * d, where m is
%   the number of vectors X. Kc is the kernel matrix of centered feature
%   vectors phi(X).
%
%   If you want to obtain the eigenvectors and eigenvalues in feature space
%   (for those cases where feature space is finite dimensional), you need
%   to compute the feature space training vectors phi(X), center them to
%   phic, and do
%
%     v = phic * a * diag(1./sqrt(da))
%     d = da / m
%
%   where m is the number of training vectors. (Note that v computed as
%   above is equal to v computed from the feature vectors save possibly the
%   sign).
%
% [1] B. Schölkopf, A.J. Smola and K.-R. Müller. "Kernel Principal
% Component Analysis". In "Advances in Kernel Methods - SV Learning", B.
% Schölkopf, C.J.C. Burges and A.J. Smola, eds. pp. 327--352, MIT Press,
% 1999.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2009-2011 University of Oxford
% Version: 1.1.0
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
error( nargchk( 1, 2, nargin, 'struct' ) );
error( nargoutchk( 0, 2, nargout, 'struct' ) );

if ndims( x ) > 2
    error( 'X must be a matrix, not a volume' )
end

% get sizes
%   N: vector length
%   M: number of training vectors
[ N, M ] = size( x );

% defaults
if ( nargin < 2 || isempty( kernel ) )
    kernel.type = 'linear';
end

switch kernel.type
    
    case 'linear'
        
        % center training vectors
        xmean = mean( x, 2 );
        x = x - repmat( xmean, 1, M );

        % compute matrix for eigenproblem
        if ( M < N ) % "kernel matrix" trick to speed up computation
            x_cov = x' * x;
        else % normal case, covariance matrix
            x_cov = x * x' / M;
        end

        % compute eigenvalues and eigenvectors
        [ v, d ] = eig( x_cov );

        % save memory by reducing eigenvalues to a vector
        d = diag( d );

        % reoder in decreasing order of signed modulus value
        [~, idx] = sort(abs(d).*sign(d), 1, 'descend');
        d = d(idx);
        v = v(:, idx);
        
        % if "kernel matrix" trick was used, it is necessary to convert the
        % coefficient eigenvectors and eigenvalues to feature space
        % eigenvectors and eigenvalues
        if ( M < N )
            
            % scale eigenvectors
            v = x * ( v * diag( 1 ./ sqrt( d ) ) );
            
            % scale eigenvalues
            d = d / M;
            
        end
        
    otherwise

        % for kernel PCA, we use a generalization of the speed up trick
        % above
        k = pts_kmat( kernel, x );
        
        % center kernel matrix (this is the same as centering the feature
        % space training vectors and then computing the kernel matrix)
        onesm = ones( M ) / M;
        k = k - onesm * k - k * onesm + onesm * k * onesm;
        
        % k must be symmetric in order to avoid complex eigenvalues, and in
        % theory it should be, but there are small errors due to numeric
        % precision when the centered kernel matrix is computed
        k = ( k + k' ) / 2;
        
        % compute eigenvalues and eigenvectors
        [ v, d ] = eig( k );

        % save memory by reducing eigenvalues to a vector
        d = diag( d );

        % reoder in decreasing order of signed modulus value
        [~, idx] = sort(abs(d).*sign(d), 1, 'descend');
        d = d(idx);
        v = v(:, idx);
        
end
