function x = elastix_transf_imcoord2(x, t)
% elastix_transf_imcoord  Convert 2D pixel coordinates using
% elastix/transformix transform.
%
% Let x be the (x, y)-coordinates of a pixel in an image (instead of
% row/column). This function computes the coordinates that the same pixel
% will have in the new image after transformix or elastix transforms.
%
% XT = elastix_transf_imcoord2(X, T)
%
%   X is a two-column matrix with the Cartesian coordinates of one or more
%   points in the input image.
%
%   T is a struct produced by elastix with the details of the image
%   transformation. See help elastix for details.
%
%   XT has the same size as X, and contains the coordinates of the same
%   points in the output image.
%
% See also: elastix, transformix, elastix_read_file2param,
% elastix_read_reg_output, elastix_write_param2file.

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

if (size(x, 2) ~= 2)
    error('X must have two columns (2D points)')
end
if (~isstruct(t))
    error('T must be a struct with the transform parameters')
end
if (isfield(t, 'Direction') && any(t.Direction ~= [1 0 0 1]))
    error('Not implemented for Direction ~= [1 0 0 1]')
end

% defaults
if (~isfield(t, 'CenterOfRotationPoint'))
    t.CenterOfRotationPoint = [0 0];
end

% select type of transform
switch (t.Transform)
    
    case 'SimilarityTransform' % similarity transform
        
        if (length(t.TransformParameters) ~= 4)
            error('SimilarityTransform must have 4 parameters')
        end
        
        % transformation nomenclature
        s =     t.TransformParameters(1);
        theta = t.TransformParameters(2);
        tx =    t.TransformParameters(3);
        ty =    t.TransformParameters(4);
        cx =    t.CenterOfRotationPoint(1);
        cy =    t.CenterOfRotationPoint(2);
        
        % transform points
        %
        % note that this is the inverse of the similarity transform. The
        % reason is the way that ITK, and by extension elastix/transformix,
        % implement image transformations
        x = [cos(theta) sin(theta); -sin(theta) cos(theta)] ...
            * [x(:, 1)' - tx - cx; x(:, 2)' - ty - cy] / s;
        x = x';
        x(:, 1) = x(:, 1) + cx;
        x(:, 2) = x(:, 2) + cy;
        
    otherwise
        
        error('Transform not implemented')
        
end
