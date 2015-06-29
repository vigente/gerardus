function t = elastix_colon(t, idx)
% ELASTIX_COLON  Colon operator for elastix transforms.
%
% This function implements the (:) operator for lists of elastix
% transforms.
%
% TOUT = ELASTIX_COLON(T, IDX)
%
%   T is a struct, or a vector of structs. Each T(i) is a list of elastix
%   transforms.
%
%   T(i) is the last transform to be applied to an image.
%   T(i).InitialTransformParametersFileName is the 2nd from last, and so
%   on. The first transform to be applied has
%   T(i).InitialTransformParametersFileName[...].InitialTransformParametersFileName   = 'NoInitialTransform'.
%
%   IDX is the list of indices. Note that index "1" means the first
%   transform, which is the last one in the list.
%
%   TOUT is T with the list of transforms modified according to IDX.
%
% Example:
%
% We have a list of transforms with 3 levels.
%
% T.Transform = 'BSplineTransform';
% T.Transform.InitialTransformParametersFileName.Transform = 'EulerTransform';
% T.Transform.InitialTransformParametersFileName.InitialTransformParametersFileName.Transform = 'TranslationTransform';
% T.Transform.InitialTransformParametersFileName.InitialTransformParametersFileName.InitialTransformParametersFileName = 'NoInitialTransform';
%
% In order to remove the very first Translation and keep the Euler and
% BSpline transforms
%
%   TOUT = ELASTIX_COLON(T, [2 3]);
%
% This produces
%
% T.Transform = 'BSplineTransform';
% T.Transform.InitialTransformParametersFileName.Transform = 'EulerTransform';
% T.Transform.InitialTransformParametersFileName.InitialTransformParametersFileName = 'NoInitialTransform';

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
narginchk(2, 2);
nargoutchk(0, 1);

% number of transforms in each level
N = length(t);

% put every transform level into a cell
tcell = {};
while (true)
    
    % extract one transform level
    tcell{end+1} = t;
    for I = 1:N
        tcell{end}(I).InitialTransformParametersFileName = 'NoInitialTransform';
    end
    
    % if this is the last level, we exit the loop
    if (strcmp(t(1).InitialTransformParametersFileName, 'NoInitialTransform'))
        break
    end
    
    % remove the top level
    t(:) = t(:).InitialTransformParametersFileName;
    
end

% invert the order of the cells so that the first transform is in the first
% cell, and the last transform is in the last cell
tcell = tcell(end:-1:1);

% keep only the indices wanted by the user
tcell = tcell(idx);

% first level of transforms
t = tcell{1};
tcell(1) = [];

while ~isempty(tcell)
    
    % add initial transform to current level
    for I = 1:N
        tcell{1}(I).InitialTransformParametersFileName = t(I);
    end
    
    % update list of transforms
    t = tcell{1};
    
    % remove current level
    tcell(1) = [];
    
end
