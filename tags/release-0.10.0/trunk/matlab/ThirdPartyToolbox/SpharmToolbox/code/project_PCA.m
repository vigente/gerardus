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

function project_PCA(fvecs, eigenvecs, eigenvals, meshsize, maxDegree, level, sigma, windowHandle, fname, cpath)

figure(windowHandle);
set(gcf,'PaperPositionMode','auto');        
% get the mean fvec
mfvec = mean(fvecs,1);


switch meshsize
    case 'quad32'
        load(fullfile(cpath,'ParamMeshes/R32_quad.mat'));
    case 'quad64'
        load(fullfile(cpath,'ParamMeshes/R64_quad.mat'));
    case 'quad128'
        load(fullfile(cpath,'ParamMeshes/R128_quad.mat'));
    case 'quad256'
        load(fullfile(cpath,'ParamMeshes/R256_quad.mat'));
    case 'quad512'
        load(fullfile(cpath,'ParamMeshes/R512_quad.mat'));
    case 'icosa1'
        load(fullfile(cpath,'ParamMeshes/L1_icosa.mat'));
    case 'icosa2'
        load(fullfile(cpath,'ParamMeshes/L2_icosa.mat'));
    case 'icosa3'
        load(fullfile(cpath,'ParamMeshes/L3_icosa.mat'));
    case 'icosa4'
        load(fullfile(cpath,'ParamMeshes/L4_icosa.mat'));
    case 'icosa5'
        load(fullfile(cpath,'ParamMeshes/L5_icosa.mat'));
    case 'icosa6'
        load(fullfile(cpath,'ParamMeshes/L6_icosa.mat'));
end

degree = sqrt(size(mfvec,2)/3)-1;
degree = min([maxDegree degree]);
Z = calculate_SPHARM_basis(sph_verts, degree);
lb = 1;
ub = (degree+1)^2;

level = min([level size(eigenvecs,2)]);

for mode = 1:level;

    eigenmode = eigenvecs(:,mode); % eigenmode direction
    variance = eigenvals(mode); % data variance in the current eigenmode
    standev = sqrt(variance);

    for i=-sigma:sigma
        % shift along the current eigenmode direction
        cfvec = mfvec + eigenmode'*i*standev;

        % visualize the current fvec
        cfvec = reshape(cfvec,3,ub)';
        spharm_vertices = real(Z(:,lb:ub)*cfvec(lb:ub,:));
        subplot(level,sigma*2+1,(mode-1)*(sigma*2+1)+i+sigma+1); 
        patch_lighta(spharm_vertices, faces); axis off;
        title(sprintf('%d*std',i));
    end
end

return;