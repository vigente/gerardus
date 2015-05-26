function tf = elastix_affine_matrix2struct(a, tf)
% ELASTIX_AFFINE_MATRIX2STRUCT  Convert homogeneous coordinates matrix
% format to elastix struct format.
%
% TF2 = ELASTIX_AFFINE_MATRIX2STRUCT(A, TF)
%
%   A is an affine transform matrix in homegeneous coordinates, referred to
%   the centre of coordinates = 0. The transform is computed as
%
%     [y 1] = [x 1] * A
%
%   where A = [a11 a12 0]
%             [a21 a22 0]
%             [t1  t2  1]
%
%   [t1 t2] is the translation vector, whereas
%
%   for a general affine transform: [a11 a12]
%                                   [a21 a22]
%
%   for a similarity transform:     [ cos(theta) sin(theta)] * s
%                                   [-sin(theta) cos(theta)]
%
%   for a rigid transform:          [ cos(theta) sin(theta)]
%                                   [-sin(theta) cos(theta)]
%
%   for a translation transform:    [1 0]
%                                   [0 1]
%
%   TF is a struct in elastix format of a 2D transform with any centre of
%   rotation. The transform may be of any affine type:
%
%     'AffineTransform'
%     'SimilarityTransform'
%     'EulerTransform'
%     'TranslationTransform'
%
%   It may also have an arbitrary TF.CenterOfRotationPoint, it doesn't need
%   to be 0.
%
%   TF2 is the output transform. It's the same as TF, but with the
%   TransformParameters recomputed, and referred to the same centre of
%   rotation.
%
% See also: elastix_affine_struct2matrix.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2015 University of Oxford
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

% split full affine matrix into blocks
t = a(3, 1:2);
A = a(1:2, 1:2);

% centre of rotation
c = tf.CenterOfRotationPoint;

% shift transform on centre of rotation
t = t - c * (eye(2) - A);

switch (tf.Transform)
    
    case 'AffineTransform'
        
        % assign parameters
        tf.TransformParameters(1) = A(1, 1);
        tf.TransformParameters(2) = A(2, 1);
        tf.TransformParameters(3) = A(1, 2);
        tf.TransformParameters(4) = A(2, 2);
        tf.TransformParameters(5:6) = t;
    
    case 'SimilarityTransform'
        
        % extract parameters
        s     = sqrt(A(1, 1).^2 + A(1, 2).^2);
        theta = atan2(A(1, 2), A(1, 1));
        
        % assign parameters
        tf.TransformParameters(1)   = s;
        tf.TransformParameters(2)   = theta;
        tf.TransformParameters(3:4) = t;
        
    case 'EulerTransform'
        
        % extract parameters
        theta = atan2(A(1, 2), A(1, 1));
        
        % assign parameters
        tf.TransformParameters(1)   = theta;
        tf.TransformParameters(2:3) = t;
        
    case 'TranslationTransform'
        
        % assign parameters
        tf.TransformParameters(1:2) = t;
        
    otherwise
        
        error('Transform not implemented')
        
end
