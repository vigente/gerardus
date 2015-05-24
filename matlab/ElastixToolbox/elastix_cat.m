function varargin = elastix_cat(varargin)
% ELASTIX_CAT  Concatenation of elastix transforms.
%
% Elastix allows for a list of transforms to be applied to an image. For
% example, if we have a translation transform
%
%   T1 = 
% 
%                              Transform: 'TranslationTransform'
%                     NumberOfParameters: 1
%                    TransformParameters: [1x2 double]
%    InitialTransformParametersFileName:  [1x1 struct]
%                 HowToCombineTransforms: 'Compose'
%
% we can see that "InitialTransformParametersFileName" points to another
% transform, and "HowToCombineTransforms" specifies that they will be
% composed. The second transform could be, for example, a B-spline
%
%   T1.InitialTransformParametersFileName = 
% 
%                              Transform: 'BSplineTransform'
%                     NumberOfParameters: 84
%                    TransformParameters: [1x84 double]
%     InitialTransformParametersFileName: 'NoInitialTransform'
%
% This second transform does not point any further transforms (although it
% could, having a longer list). In this example, the translation is applied
% first to the image, and then the B-spline.
%
% Note: It's a bit confusing that "InitialTransformParametersFileName"
% points to the transform that is going to be applied to the image after
% the current one. However, this makes sense because when we register A to
% B in elastix, the transform returned is for coordinates from B to A. The
% reason for this is that otherwise transforming an image could leave
% "holes" after resampling.
%
% ELASTIX_CAT allows to concatenate a list of elastix transforms, each of
% which can be a simple transform or a list of transforms.
%
% TOUT = ELASTIX_CAT(T1, ..., TN)
%
%   T1, ..., TN is a list of N elastix transforms, each being a struct.
%   Each transform can be simple, or be a concatenation of transforms using
%   the InitialTransformParametersFileName field.
%
%   TOUT is a struct with T1, ..., TN concatenated such that T1 is applied
%   first to the image, then T2, etc, until TN.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2015 University of Oxford
% Version: 0.1.1
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
narginchk(1, Inf);
nargoutchk(0, 1);

% number of transforms to concatenate
N = length(varargin);

if (N == 0)
    
    varargin = [];
    return;
    
elseif (N == 1)
    
    varargin = varargin{1};
    
end

% loop from the back of the list of transforms
for I = N-1:-1:1
    
    varargin{I} = cat_2_transf(varargin{I}, varargin{I+1});
    
end
varargin = varargin{1};

end

%% local functions

% cat_2_transf: add transform t2 at the end of transform t1
function t1 = cat_2_transf(t1, t2)

% look for the end 
str = 't1';
while (isstruct(eval([str '.InitialTransformParametersFileName'])))
    str = [str '.InitialTransformParametersFileName'];
end
eval([str '.InitialTransformParametersFileName = t2;']);

end
