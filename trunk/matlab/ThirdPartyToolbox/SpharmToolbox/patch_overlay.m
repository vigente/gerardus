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

function patch_overlay(vertices, faces, signal, threshold, type, vtnorm, cmap,windowHandle, fname)

figure(windowHandle);
set(gcf,'PaperPositionMode','auto');
colormap(cmap);

patch_lighta(vertices, faces);
maxSG = max(signal);
minSG = min(signal);

FValpha = ones(size(vertices,1),1);
switch type
    case 'p-value'
        idx = find(abs(signal)>(threshold));
        FValpha(idx) = 0;

        hold on; 
        patch('faces', faces,'vertices',vertices, 'FaceVertexCData',signal, 'EdgeColor', 'none', ... 
            'VertexNormals',vtnorm,'FaceVertexAlpha',FValpha,'FaceColor','interp', 'FaceAlpha','interp'); 
        hold off;
        colorbar;
        
    case 't-map'
        idx = find(abs(signal)<threshold);
        FValpha(idx) = 0;

        hold on; 
        patch('faces', faces,'vertices',vertices, 'FaceVertexCData',signal, 'EdgeColor', 'none', ... 
            'VertexNormals',vtnorm,'FaceVertexAlpha',FValpha,'FaceColor','interp', 'FaceAlpha','interp'); 
        hold off;
        colorbar('location', 'EastOutside');
end

title(sprintf('%s.mat\n%s',strrep(fname,'_','-'),type));

return;