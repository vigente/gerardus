function h = fspecial3(type, sz)
% FSPECIAL3  Create predefined 3-dimensional filters
%
% H = FSPECIAL3(TYPE, SZ)
%
%   H is a 3-dimensional (3D) filter of TYPE:
%
%     'gaussian'  Rotationally symmetric Gaussian low-pass filter (default)
%
%   SZ is a vector with the size of the output filter in each dimension, in
%   order [rows, columns, slices]. By default, SZ = [3 3 3]. Sizes have to
%   be odd numbers so that the filter can be centered around 0.
%
% See also: fspecial.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
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
error(nargchk(0, 2, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

% defaults
if (nargin < 1 || isempty(type))
    type = 'gaussian';
end
if (nargin < 2 || isempty(sz))
    sz = [3 3 3];
end

% check that sizes are odd numbers
if (any(~rem(sz,2)))
    error('Sizes have to be odd numbers')
end

% filter size to each side of 0
l = floor(sz/2);

% compute filter
switch type
    
    case 'gaussian'
        
        % filter domain
        [gr, gc, gs] = ndgrid(-l(1):l(1), -l(2):l(2), -l(3):l(3));
        
        % compute gaussian function
        h = exp(-(gr.*gr + gc.*gc + gs.*gs)*.5);
        
        % normalize filtered intensity
        h = h / sum(h(:));

    otherwise
        
        error('Filter type not implemented')
end
