% test_ventricles_phantom.m
%
% Script to test the generation of a heart phantom

% Author: Christopher Kelly <christopher.kelly28@gmail.com>
% Copyright © 2013 University of Oxford
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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

%% run function with zero noise

[heart,LV,RV,myo] = ventricles_phantom;

%% view result with no noise

x_range = -(size(heart,1)/2)+0.5:(size(heart,1)/2)-0.5;
y_range = -(size(heart,2)/2)+0.5:(size(heart,2)/2)-0.5;
z_range = -(size(heart,3)/2)+0.5:(size(heart,3)/2)-0.5;
[x,y,z] = meshgrid(x_range,y_range,z_range);

figure;view(3)
surface( x(:,:,size(heart,3)/2),...
         y(:,:,size(heart,3)/2),...
         z(:,:,size(heart,3)/2),...
         heart(:,:,size(heart,3)/2),...
         'edgecolor', 'none');

surface( permute(x(:,size(heart,1)/2,:),[1,3,2]),...
         permute(y(:,size(heart,1)/2,:),[1,3,2]),...
         permute(z(:,size(heart,1)/2,:),[1,3,2]),...
         permute(heart(:,size(heart,1)/2,:),[1,3,2]),...
         'edgecolor', 'none'); 
     
surface( permute(x((size(heart,2)/2),:,:),[2,3,1]),...
         permute(y((size(heart,2)/2),:,:),[2,3,1]),...
         permute(z((size(heart,2)/2),:,:),[2,3,1]),...
         permute(heart((size(heart,2)/2),:,:),[2,3,1]),...
         'edgecolor', 'none');
     
%% run function with noise added

[heart,LV,RV,myo] = ventricles_phantom(0.05,0.03);

%% view result with noise

figure;view(3)
surface( x(:,:,size(heart,3)/2),...
         y(:,:,size(heart,3)/2),...
         z(:,:,size(heart,3)/2),...
         heart(:,:,size(heart,3)/2),...
         'edgecolor', 'none');

surface( permute(x(:,size(heart,1)/2,:),[1,3,2]),...
         permute(y(:,size(heart,1)/2,:),[1,3,2]),...
         permute(z(:,size(heart,1)/2,:),[1,3,2]),...
         permute(heart(:,size(heart,1)/2,:),[1,3,2]),...
         'edgecolor', 'none'); 
     
surface( permute(x((size(heart,2)/2),:,:),[2,3,1]),...
         permute(y((size(heart,2)/2),:,:),[2,3,1]),...
         permute(z((size(heart,2)/2),:,:),[2,3,1]),...
         permute(heart((size(heart,2)/2),:,:),[2,3,1]),...
         'edgecolor', 'none');
     
%% create a mesh from the myocardium binary mask

% run function with zero noise
[heart, LV, RV, myo] = ventricles_phantom;

% convert myocardial binary mask to mesh at full resolution
[x, tri] = binsurface(myo, 3);

% plot mesh
hold off
plotmesh(x, tri)

% convert myocardial binary mask to mesh at reduced resolution
isovalues = 1;
% opt = 1.5;
% method = 'cgalsurf';
opt = .1;
method = 'simplify';
tic
[x, tri] = v2s(single(myo), isovalues, opt, method);
tri = tri(:, 1:3);
toc

% plot mesh
cla
plotmesh(x, tri)
view(-62, 52)
axis equal
