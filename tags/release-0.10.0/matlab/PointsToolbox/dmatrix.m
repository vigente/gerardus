function d = dmatrix(x, y, METHOD)
% DMATRIX  Matrix of distances between vectors
%
% D = DMATRIX(X)
%
%   D(i,j) is the Euclidean distance between X(:,i) and X(:,j)
%
% D = DMATRIX(X, Y)
%
%   D(i,j) is the Euclidean distance between X(:,i) and Y(:,j)
%
% D = DMATRIX(X, Y, METHOD)
%
%   METHOD is a string with the type of distance measured. By default,
%   METHOD='euclidean'. METHOD='centrip' computes the centripetal
%   distance [1].
%
% [1] ETY Lee, "Choosing nodes in parametric curve interpolation",
% Computer-Aided Design, 21:363–370, 1989.

% Author: Ramon Casero <rcasero@gmail.com>.
% Copyright © 2011 University of Oxford
% v0.3.1
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
narginchk(1, 3);
nargoutchk(0, 1);

% defaults
if (nargin < 3 || isempty(METHOD))
    METHOD='euclidean';
end

%    D = DMATRIX(X)
if (nargin == 1)
    
    M = size(x, 2);
    d = zeros(M);
    for I = 1:M
        for J = (I + 1):M
            d(I, J) = sqrt(sum((x(:, I) - x(:, J)) .^ 2));
            d(J, I) = d(I, J);
        end
    end

%    D = DMATRIX(X,Y)
else
    
    M = size(x, 2);
    N = size(y, 2);
    d = zeros(M, N);
    for I = 1:M
        for J = 1:N
            d(I, J) = sqrt(sum((x(:, I) - y(:, J)) .^ 2));
        end
    end

end

% centripetal distance instead of Euclidean?
if strcmp(METHOD, 'centrip')
    d = sqrt(d);
end
