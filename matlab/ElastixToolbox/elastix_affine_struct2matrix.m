function a = elastix_affine_struct2matrix(tf)
% ELASTIX_AFFINE_STRUCT2MATRIX  Convert any type of affine transform from
% elastix struct format to homogeneous coordinates matrix format.
%
% A = ELASTIX_AFFINE_STRUCT2MATRIX(TF)
%
%   TF is a struct in elastix format of a 2D transform. The transform may
%   be of any affine type:
%
%     'AffineTransform'
%     'SimilarityTransform'
%     'EulerTransform'
%     'TranslationTransform'
%
%   It may also have an arbitrary TF.CenterOfRotationPoint, it doesn't need
%   to be 0. Subsequent transforms pointed at from
%   TF.InitialTransformParametersFileName are ignored.
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
% See also: elastix_affine_matrix2struct.

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
narginchk(1, 1);
nargoutchk(0, 1);

switch (tf.Transform)
    
    case 'AffineTransform'
        
        % transformation nomenclature
        a11 =   tf.TransformParameters(1);
        a21 =   tf.TransformParameters(2);
        a12 =   tf.TransformParameters(3);
        a22 =   tf.TransformParameters(4);
        t =     tf.TransformParameters(5:6);
    
        % matrix block that doesn't contain the translation
        A = [a11 a12; a21 a22];
        
    case 'SimilarityTransform'
        
        % transformation nomenclature
        s =     tf.TransformParameters(1);
        theta = tf.TransformParameters(2);
        t =     tf.TransformParameters(3:4);
        
        % matrix block that doesn't contain the translation
        A = [cos(theta) sin(theta);...
            -sin(theta) cos(theta)] * s;
        
    case 'EulerTransform'
        
        % transformation nomenclature
        theta = tf.TransformParameters(1);
        t =     tf.TransformParameters(2:3);
        
        % matrix block that doesn't contain the translation
        A = [cos(theta) sin(theta);...
            -sin(theta) cos(theta)];
        
    case 'TranslationTransform'
        
        % transformation nomenclature
        t =     tf.TransformParameters(1:2);
        
        % matrix block that doesn't contain the translation
        A = [1 0;...
            0 1];
        
    otherwise
        
        error('Transform not implemented')
        
end

% centre of rotation
c = tf.CenterOfRotationPoint;

% center transform on origin
t = t + c * (eye(2) - A);

% create affine transform matrix
a = [A [0;0]; t 1];
