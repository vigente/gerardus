% test_elastix_compose_transf.m

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Combination of similarity transforms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% toy image 
im = zeros(50, 76, 'uint8');
im(1:25, 1:38) = 1;

% load example transform
load('SimilarityTransformExample.mat', 't')

% change default pixel value, so that we see where new background is being
% created
t.DefaultPixelValue = 2;
t.FinalBSplineInterpolationOrder = 0;

% create two transforms
t1 = t;
t2 = t;

% plot output
subplot(2, 2, 1)
imagesc(im)
axis xy equal

% first transform
t1.CenterOfRotationPoint = [38 25];
t1.TransformParameters = [2 30/180*pi 20 10];

% transform image
im2 = transformix(t1, im);

% plot output
subplot(2, 2, 2)
imagesc(im2)
axis xy equal

% second transform
t2.CenterOfRotationPoint = [30 40];
t2.TransformParameters = [1.5 80/180*pi -15 -10];

% transform im2
im3 = transformix(t2, im2);

% plot output
subplot(2, 2, 3)
imagesc(im3)
axis xy equal

% compute combined transform
tc = elastix_compose_transf(t1, t2);

% transform im with the combined transform
imc = transformix(tc, im);

% plot output
subplot(2, 2, 4)
imagesc(imc)
axis xy equal

% apply the two transforms to the point on the top right corner
aux = elastix_transf_imcoord2([37 24], t1);
elastix_transf_imcoord2(aux, t2)

% apply composed transform to the same point (both results should be the
% same)
elastix_transf_imcoord2([37 24], tc)
