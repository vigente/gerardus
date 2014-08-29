function tfc = elastix_compose_transf(tf1, tf2)
% elastix_compose_transf  Composition of two 2D affine transforms.
%
% elastix_compose_transf composes two 2D transforms from the affine family,
% and produces another affine transform.
%
% TFC = elastix_compose_transf(TF1, TF2)
%
%   TF1, TF2 are two transforms, either given as (3, 3)-matrices with the
%   Matlab tform convention (help affine2d and projective2d for details) to
%   map coordinates in homogeneous coordinates with the rotation centered
%   on the origin of coordinates (0, 0):
%
%   [y 1] = [x 1] * [A 0]
%                   [t 1]
%
%   Note that for rigid transformations, A = [ cos(theta)  sin(theta)]
%                                            [-sin(theta)  cos(theta)]
%
%   or in the struct format produced by elastix (see help elastix for
%   details):
%
%   'AffineTransform'
%   'SimilarityTransform (Currently only implemented for this one).
%   'EulerTransform' (Rigid)
%   'TranslationTransform'
%
%   which can be centered anywhere else.
%
%   TFC is the composed transform. TFC is a matrix or struct depending on
%   whether TFC1 is a matrix or struct, and it has the same centre of
%   rotation. If TF1 is a struct, then TFC is given as the most flexible of
%   the two types in TF1, TF2. For example, if TF1 is EulerTransform and
%   TF2 is SimilarityTransform, TFC will be SimilarityTransform.
%
%   transformix(TFC, IM) produces the same result in one step than
%   transformix(TF2, transformix(TF1, IM)).
%
%
% See also: elastix, transformix, elastix_transf_imcoord2.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
% Version: 0.2.1
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

if (~strcmp(tf1.HowToCombineTransforms, 'Compose') ...
        || ~strcmp(tf2.HowToCombineTransforms, 'Compose'))
    error('HowToCombineTransforms must be ''Compose''')
end

% convert transforms to affine matrix with center of rotation = 0
a1 = origin_affine(tf1);
a2 = origin_affine(tf2);

% compose transforms
ac = a2 * a1;

% format output transform and center on same center as first transform
tfc = format_centered_output(ac, tf1, tf2);

end

function a = origin_affine(tf)

% transform is provided as an elastix struct
if (isstruct(tf))
    
    switch (tf.Transform)
        
        case 'SimilarityTransform'
            
            % transformation nomenclature
            s =     tf.TransformParameters(1);
            theta = tf.TransformParameters(2);
            t =     tf.TransformParameters(3:4);
            c =     tf.CenterOfRotationPoint;
            
            % center transform on origin
            R = [cos(theta) sin(theta);...
                -sin(theta) cos(theta)];
            t = t + c * (eye(2) - s * R);
            
            % create affine transform matrix
            a = [s*R [0;0]; t 1];
            
        otherwise
            
            error('Transform not implemented')
            
    end
    
else % transform is provided as an affine matrix
    
    if (any(size(tf) ~= [3 3]))
        error('If TF is provided as an affine matrix, it must be a (3, 3)-matrix')
    end
    
    % we don't need to do anything to the matrix
    a = tf;
    
end

end

function tfc = format_centered_output(ac, tf1, tf2)

if (isstruct(tf1)) % return composed transform as an elastix struct
    
    % use TF1 as a template for the output
    tfc = tf1;
    
    % the return type is the most flexible of the two transforms
    if (strcmp(tf1.Transform, 'AffineTransform') ...
            || strcmp(tf2.Transform, 'AffineTransform'))
        
        error('AffineTransform not implemented')
        
    elseif (strcmp(tf1.Transform, 'SimilarityTransform') ...
            || strcmp(tf2.Transform, 'SimilarityTransform'))

        % convert general affine matrix into similarity transform
        % parameters
        s = sqrt(ac(1, 1)^2 + ac(1, 2)^2);
        t = ac(3, 1:2);
        theta = atan2(ac(1, 2), ac(1, 1));
        
        % center transform on same center as TF1
        c = tf1.CenterOfRotationPoint;
        R = [ cos(theta) sin(theta);...
             -sin(theta) cos(theta)];
        t = t - c * (eye(2) - s * R);

        % assign parameters to output struct
        tfc.TransformParameters(1) = s;
        tfc.TransformParameters(2) = theta;
        tfc.TransformParameters(3:4) = t;
        tfc.CenterOfRotationPoint = c;
        
    elseif (strcmp(tf1.Transform, 'EulerTransform') ...
            || strcmp(tf2.Transform, 'EulerTransform'))
        
        error('EulerTransform (Rigid) not implemented')
        
    elseif (strcmp(tf1.Transform, 'TranslationTransform') ...
            || strcmp(tf2.Transform, 'TranslationTransform'))
        
        error('TranslationTransform not implemented')
        
    else
            
            error('Unknow transform in TF1 or TF2')
            
    end
    
else % return composed transform as a matrix
    
    tfc = ac;
    
end

end
