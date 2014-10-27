function v = elastix_affine_tform2param(tform, type)
% elastix_affine_tform2param  Convert Matlab tform affine transform to
% elastix parameter vector.
%
% V = elastix_affine_tform2param(TFORM, TYPE)
%
%   TFORM is a Matlab tform with a TYPE='translation', 'rigid',
%   'similarity' or 'affine' transform.
%
%   Currently, only similarity transform is implemented.
%
%   V is a vector with the transform parameters in the format used by
%   elastix and ITK:
%
%     'similarity': [s theta tx ty]
%
%   If TFORM is a vector of structs, then V is a matrix, with one row per
%   element in TFORM.

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
narginchk(2, 2);
nargoutchk(0, 1);

switch type
    
    case 'similarity'
        
        % init output
        v = zeros(length(tform), 4);
        
        for I = 1:length(tform)

            % extract transform parameters
            s = sqrt(sum(tform(I).T(1, 1:2) .^2));
            theta = atan2(tform(I).T(1, 2), tform(I).T(1, 1));
            tx = tform(I).T(3, 1);
            ty = tform(I).T(3, 2);
            
            % format output
            v(I, :) = [s theta tx ty];
            
        end
        
    otherwise
        
        error('TYPE not implemeted')
        
end
