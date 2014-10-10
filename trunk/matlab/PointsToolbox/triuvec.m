function x = triuvec(x, k)
% triuvec  Extract upper triangular part into vector form.
%
% Matlab's triu() returns a matrix, with all elements below the upper
% triangular part equal to 0. triuvec() extracts only the elements in the
% upper triangular part, and returns them as a vector.
%
% Y = triuvec(X)
%
%   X is a matrix.
%
%   Y is the upper triangular part of X in vector form.
%
% Y = triuvec(X, K)
%
%   K selects another diagonal. K>0 is above the main diagonal and K<0 is
%   below the main diagonal.
%
% See also: triu.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
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
narginchk(1, 2);
nargoutchk(0, 1);

% defaults
if (nargin < 2)
    k = 0;
end

% extract upper triangular part into vector form
x = x(triu(true(size(x)), k));
