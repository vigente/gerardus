function tfc = elastix_compose_afftransf(tf1, tf2)
% ELASTIX_COMPOSE_AFFTRANSF  Composition of two 2D affine transforms.
%
% ELASTIX_COMPOSE_AFFTRANSF composes two 2D transforms from the affine
% family, and produces another affine transform.
%
% TFC = ELASTIX_COMPOSE_AFFTRANSF(TF1, TF2)
%
%   TF1, TF2 are two affine transforms to be applied in that order to an
%   _image_ (not to point coordinates). E.g. TF1 is a translation
%   transform, and TF2 is a rigid transform.
%
%   IMPORTANT: In elastix,
%
%       TF2.InitialTransformParametersFileName = TF1
%
%   means that TF1 is the initial transform, followed by TF2. However, if
%   you want to manually apply the transforms to the image, you have to do
%   in the inverse order to obtain the same result (except for accumulated
%   interpolation errors):
%
%       IM2 = transformix(TF2, IM);
%       IM3 = transformix(TF1, IM2);
%
%   TF1, TF2 can have the following formats:
%
%   * (3, 3)-matrices with the Matlab tform convention (help affine2d and
%     projective2d for details) to map coordinates in homogeneous
%     coordinates with the rotation centered on the origin of coordinates
%     (0, 0):
%
%     [y 1] = [x 1] * [A 0]
%                     [t 1]
%
%     Note that for rigid transformations, A = [ cos(theta)  sin(theta)]
%                                              [-sin(theta)  cos(theta)]
%
%   * struct format produced by elastix (see help elastix for details):
%
%     'AffineTransform'
%     'SimilarityTransform'
%     'EulerTransform' (= Rigid transform)
%     'TranslationTransform'
%
%     which can be centered anywhere else. Note that if TF1 or TF2 point to
%     any subsequent transforms with TF.InitialTransformParametersFileName,
%     they will be ignored.
%
%   TFC is the composed transform. TFC is a matrix or struct depending on
%   whether TFC1 is a matrix or struct, and it has the same centre of
%   rotation as TF1. If TF1 is a struct, then TFC is given as the most
%   flexible of the two types in TF1, TF2. For example, if TF1 is
%   an EulerTransform and TF2 is a SimilarityTransform, TFC will be
%   a SimilarityTransform.
%
%   TFC has the same output size as TF2, which is logical, because TF2 is
%   the last transform to be applied to the image. However, we keep the
%   center of rotation of TF1.
%
%
% See also: elastix, transformix, elastix_transf_imcoord2.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014-2015 University of Oxford
% Version: 0.2.6
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

% if one of the transforms is empty, the composition returns the other
% transform
if (isempty(tf1) || (isstruct(tf1) && isempty(tf1.Transform)))
    tfc = tf2;
    return
elseif (isempty(tf2) || (isstruct(tf2) && isempty(tf2.Transform)))
    tfc = tf1;
    return
end

if ((isstruct(tf1) && ~strcmp(tf1.HowToCombineTransforms, 'Compose')) ...
        || (isstruct(tf2) && ~strcmp(tf2.HowToCombineTransforms, 'Compose')))
    error('HowToCombineTransforms must be ''Compose''')
end

if (isstruct(tf1.InitialTransformParametersFileName))
    warning('TF1 is pointing to another transform in InitialTransformParametersFileName, that will be ignored')
    tf1.InitialTransformParametersFileName = 'NoInitialTransform';
end
if (isstruct(tf2.InitialTransformParametersFileName))
    warning('TF2 is pointing to another transform in InitialTransformParametersFileName, that will be ignored')
    tf2.InitialTransformParametersFileName = 'NoInitialTransform';
end

% convert transforms to affine matrix with center of rotation = 0
a1 = elastix_affine_struct2matrix(tf1);
a2 = elastix_affine_struct2matrix(tf2);

% compose transforms
ac = a1 * a2;

% format output transform and center on same center as first transform
if (isstruct(tf1)) % return composed transform as an elastix struct
    
    % use TF2 as a template for the output, because TF2 is the last
    % transform to be applied to the image
    tfc = tf2;
    
    % but we use the rotation center of the first transform
    tfc.CenterOfRotationPoint = tf1.CenterOfRotationPoint;
    
    % the return type is the most flexible of the two transforms
    if (strcmp(tf1.Transform, 'AffineTransform') ...
            || strcmp(tf2.Transform, 'AffineTransform'))
        
        tfc.Transform = 'AffineTransform';
        
    elseif (strcmp(tf1.Transform, 'SimilarityTransform') ...
            || strcmp(tf2.Transform, 'SimilarityTransform'))

        
        tfc.Transform = 'SimilarityTransform';
        
    elseif (strcmp(tf1.Transform, 'EulerTransform') ...
            || strcmp(tf2.Transform, 'EulerTransform')) % rigid transform
        
        tfc.Transform = 'EulerTransform';
        
    elseif (strcmp(tf1.Transform, 'TranslationTransform') ...
            || strcmp(tf2.Transform, 'TranslationTransform'))
        
        tfc.Transform = 'TranslationTransform';
        
    else
            
        error('Unknow transform in TF1 or TF2')
            
    end

    % convert full affine matrix to elastix format
    tfc = elastix_affine_matrix2struct(ac, tfc);

else % return composed transform as a matrix
    
    tfc = ac;
    
end
