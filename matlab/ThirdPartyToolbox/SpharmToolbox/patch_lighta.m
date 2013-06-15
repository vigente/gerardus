%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherical Harmonic Modeling and Analysis Toolkit (SPHARM-MAT) is a 3D 
% shape modeling and analysis toolkit. 
% It is a software package developed at Shenlab in Center for Neuroimaging, 
% Indiana University (SpharmMat@gmail.com, http://www.iupui.edu/~shenlab/)
% It is available to the scientific community as copyright freeware 
% under the terms of the GNU General Public Licence.
% 
% Copyright 2009, 2010, ShenLab, Center for Neuroimaging, Indiana University
% 
% This file is part of SPHARM-MAT.
% 
% SPHARM-MAT is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% SPHARM-MAT is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with SPHARM-MAT. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ============================================
% patch_light.m
%
% Goal: surface rendering with lighting
%       need to open figure or subplot before use
%
% Li Shen 
% 02/21/2006 - create

function patch_light(vertices,faces)

patch('faces', faces, 'vertices', vertices, ...
	'FaceVertexCData', ones(size(vertices,1),1)*[.8 .8 .8], ...
	'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 1, 'AmbientStrength',0.8); % transparancy can be adjusted by FaceAlpha
l1 = light('position', [-1 -1 .1], 'color', [.3 .3 .3]*2);
l2 = light('position', [1 -1 .1], 'color', [.3 .3 .3]*2);
l3 = light('position', [-.5 1 .1], 'color', [.3 .3 .3]*2);
l4 = light;
l5 = light('Position',[-1 0 0],'Style','infinite');
% material shiny;
% material metal;
material([.3 .4 .2 10]);
%  	lighting phong;
lighting gouraud	
axis image; box on; view(3);
    
return;

%
% surface mesh rendering
%

function patch_mesh(vertices,faces)

patch('faces', faces, 'vertices', vertices, ...
	'FaceColor', 'w', 'EdgeColor', 'k', 'FaceAlpha', 1); % transparancy can be adjusted by FaceAlpha
axis image; box on; view(3);
    
return
