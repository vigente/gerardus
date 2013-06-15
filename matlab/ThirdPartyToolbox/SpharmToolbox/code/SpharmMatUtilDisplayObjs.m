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

function SpharmMatUtilDisplayObjs(confs, objs, cpath)

numSbj = size(objs,2);

tplt = [];
if exist(confs.Template,'file')
    if strcmp(confs.Template(end-3:end),'.mat')==1
        tplt = load(confs.Template);
        if ~isfield(tplt,'fvec')
            disp('Template file does not contain fvec');
            tplt = [];
        end
    else
        disp('Template should be a .mat file');
    end
end

for i = 1:numSbj
    file = objs{i};
    [pa,na,ex,ve]=fileparts(file);
    nastr = strrep(na,'_','-');
    
    clear adc; clear cdata; rmsd = [];
    
    suffix = '';

    load(file);
    if exist('bim','var')
        roi = bim;
        DIM = size(roi);

        % Make binary image
        ind = find(roi>0); roi(ind) = 1;
        ind = find(roi<1); roi(ind) = 0;    

        % Generate voxel surface for the binary image
        [vertices, faces] =  gen_surf_data(roi,origin,vxsize);
    end
        
    h = figure('Name', ['Object: ' na],'NumberTitle', 'off');
    cameratoolbar(h, 'Show');
    
    % calculate SPHARM reconstruction
    if exist('fvec','var')
        % get regular mesh in the parameter space: sph_verts, faces
        switch strtrim(char(confs.Mesh))
            % to-do-item: add an option for reconstruction using original mesh
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
        

        
        % calculate SPHARM reconstruction using the regular mesh: vertices
        if ~strcmp(deblank(char(confs.Mesh)), 'orig')            
            % get the user-specified degree
            if isempty(confs.Degree)
                degree = Inf;
            else
                degree = confs.Degree;
            end           
            degree = min(degree,sqrt(size(fvec,1))-1);
            suffix = sprintf('%s_d%d',strtrim(char(confs.Mesh)),degree);
            % do the calculation
            Z = calculate_SPHARM_basis(sph_verts, degree);
            lb = 1;
            ub = (degree+1)^2;
            vertices = real(Z(:,lb:ub)*fvec(lb:ub,:));
            % calculate RMSD if there is a template
            if ~isempty(tplt)
                tfvec = tplt.fvec;
                tcnt = size(tfvec,1);
                if tcnt>=ub
                    tfvec = tfvec(1:ub,:);
                else
                    tfvec = zeros(1:tcnt,3);
                    tfvec(1:tcnt,:) = tplt.fvec;
                end
                rmsd = SPHARM_rmsd(fvec,tfvec);
            end
        end
    end
    
    % option for display in object space, parameter space, or both
    switch strtrim(char(confs.Space))
        case 'object'
            if ~exist('vertices','var') | ~exist('faces','var')
                disp(sprintf('Cannot display %s in the object space',nastr));
                break;
            end
            subp = 0;
            titleStr = sprintf('Object: %s',nastr);
        case 'param'
            if ~exist('sph_verts','var') | ~exist('faces','var')
                disp(sprintf('Cannot display %s in the parameter space',nastr));
                break;
            end
            subp = 1;
            titleStr = sprintf('Parameterization: %s',nastr);
        case 'both'
            if ~exist('faces','var') | (~exist('vertices','var') & ~exist('sph_verts','var'))
                disp(sprintf('Cannot display %s in object and parameter spaces',nastr));
                break;
            end
            subp = 2;
            titleStrLeft = sprintf('%s',nastr);
            titleStrRight = sprintf('Parameterization');           
    end
    
    % option to visualize parameterization in color
    % ADC is shown together with this option
    over1 = 0;
    switch strtrim(char(confs.Overlay))
        case 'none'
            overl = 0;
        case 'adc_paramap'
            if exist('vertices','var') & exist('sph_verts','var') & exist('faces','var')
                overl = 1;
                % calculate the color coding
                [THETA,PHI,R] = cart2sph(sph_verts(:,1),sph_verts(:,2),sph_verts(:,3));
                ix = find(THETA<0); THETA(ix) = THETA(ix)+pi*2;
                cdata = THETA.*sign(PHI);
                % calculate the area distortion cost: adc
                asr = calc_asr(vertices, sph_verts, faces); adc = asr.prmstc;
            else
                disp(sprintf('Color coded parameterization cannot be calculated for %s due to missing variables',nastr));
            end
    end

    % option for solid shading, mesh view, or both
    switch strtrim(char(confs.Shade))
        case 'solid'
            if overl == 1
                if subp == 0
                    patch_color(vertices, faces, cdata); title(titleStr);
                elseif subp == 1
                    patch_color(sph_verts, faces, cdata); title(titleStr);
                elseif subp == 2
                    subplot(1,2,1);
                    patch_color(vertices, faces, cdata); title(titleStrLeft);
                    subplot(1,2,2);
                    patch_color(sph_verts, faces, cdata); title(sprintf('%s: ADC=%f',titleStrRight,adc));
                end
            else
                if subp == 0
                    patch_light(vertices, faces); title(titleStr);
                elseif subp == 1
                    patch_light(sph_verts, faces); title(titleStr);
                elseif subp == 2
                    subplot(1,2,1);
                    patch_light(vertices, faces); title(titleStrLeft);
                    subplot(1,2,2);
                    patch_light(sph_verts, faces); title(titleStrRight);
                end
            end
        case 'mesh'
            if overl == 1
                disp('Color coded parameterization is not shown in mesh view');
                if subp == 0
                    patch_mesh(vertices, faces); title(titleStr);
                elseif subp == 1
                    patch_mesh(sph_verts, faces); title(titleStr);
                elseif subp == 2
                    subplot(1,2,1);
                    patch_mesh(vertices, faces); title(titleStrLeft);
                    subplot(1,2,2);
                    patch_mesh(sph_verts, faces); title(sprintf('%s: ADC=%f',titleStrRight,adc));
                end
            else
                if subp == 0
                    patch_mesh(vertices, faces); title(titleStr);
                elseif subp == 1
                    patch_mesh(sph_verts, faces); title(titleStr);
                elseif subp == 2
                    subplot(1,2,1);
                    patch_mesh(vertices, faces); title(titleStrLeft);
                    subplot(1,2,2);
                    patch_mesh(sph_verts, faces); title(titleStrRight);
                end
            end
        case 'both'
            if overl == 1
                if subp == 0
                    patch_colormesh(vertices, faces, cdata); title(titleStr);
                elseif subp == 1
                    patch_colormesh(sph_verts, faces, cdata); title(titleStr);
                elseif subp == 2
                    subplot(1,2,1);
                    patch_colormesh(vertices, faces, cdata); title(titleStrLeft);
                    subplot(1,2,2);
                    patch_colormesh(sph_verts, faces, cdata); title(sprintf('%s: ADC=%f',titleStrRight,adc));
                end
            else
                if subp == 0
                    patch_lightmesh(vertices, faces); title(titleStr);
                elseif subp == 1
                    patch_lightmesh(sph_verts, faces); title(titleStr);
                elseif subp == 2
                    subplot(1,2,1);
                    patch_lightmesh(vertices, faces); title(titleStrLeft);
                    subplot(1,2,2);
                    patch_lightmesh(sph_verts, faces); title(titleStrRight);
                end
            end
    end
    
    if ~isempty(suffix)
        suffix = sprintf('%s_%s',confs.Space(1:3),suffix); 
    else
        suffix = confs.Space(1:3);
    end
    if ~isempty(rmsd)
        zstr = sprintf('%s, rmsd=%0.3f',suffix,rmsd);
        zlabel(strrep(zstr,'_','-'));
    else
        zlabel(strrep(suffix,'_','-'));
    end
    
    % display and/or export
    switch strtrim(char(confs.Export))
        case 'screen'
            % display the figure 
        case 'png'
            % save as png 
            pa_png = fullfile(pa,'PNG');
            if ~exist(pa_png,'dir')
                mkdir(pa_png);
            end
            saveas(h,fullfile(pa_png,[na '_' suffix]),'png'); close(h);
        case 'both'
            % do both
            pa_png = fullfile(pa,'PNG');
            if ~exist(pa_png,'dir')
                mkdir(pa_png);
            end
            saveas(h,fullfile(pa_png,[na '_' suffix]),'png');            
    end
        
end

return;
