function d2 = gram2dsq(g)
% GRAM2DSQ  Convert Gram matrix to squared distance matrix
%
% D2 = GRAM2DSQ(G)
%
%   G is a Gram matrix. The Gram matrix of matrix Y, where each column has
%   the coordinates of a point, is G = Y'*Y.
%
%   D2 is the matrix of squared Euclidean distances between the points in
%   Y. Note that if Y is available, D2 = DMATRIX(Y).^2 = GRAM2DSQ(Y'*Y).

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.


% check arguments
error(nargchk(1, 1, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% compute diagonal of Gram matrix
delta = diag(g);

% compute auxiliary step
aux = delta * ones(1, size(g, 1));

% compute squared distance matrix
d2 = aux + aux' - 2*g;
