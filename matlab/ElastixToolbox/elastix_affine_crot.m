function tElx = elastix_affine_crot(tElx, c2)
% ELASTIX_AFFINE_CROT  Change centre of rotation of Elastix affine
% transform.
%
% TOUT = ELASTIX_AFFINE_CROT(T, C)
%
%   This function changes the centre of rotation of affine transform T to
%   vector C, producing TOUT. Both T and TOUT produce the same transform on
%   an image. By default, C=0.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2015 University of Oxford
% Version: 0.1.0
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
    c2 = zeros(size(tElx.CenterOfRotationPoint));
end

% compute affine matrix in homogeneous coordinates
a = elastix_affine_struct2matrix(tElx);

% nomenclature
A = a(1:end-1, 1:end-1);
t = a(end, 1:end-1);
c1 = tElx.CenterOfRotationPoint;

% update translation vector
N = size(a, 1) - 1;
t = t + (c1 - c2) * (eye(N) - A);

% update affine matrix
a(end, 1:end-1) = t;

% convert back to elastix format
tElx = elastix_affine_matrix2struct(a, tElx);

% change centre of rotation
tElx.CenterOfRotationPoint = c2;
