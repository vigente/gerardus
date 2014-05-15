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

function SpharmMatUtilFormatConvert(confs, objs, method)

numSbj = length(objs);
if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end
if (~exist([confs.OutDirectory '/Logs'],'dir'))
    mkdir([confs.OutDirectory '/Logs']);
end

%class(confs) - struct
%confs.vars -list of variable
%class(objs) - cell

h = waitbar(0,'Please wait...');
for i = 1:numSbj
    file = objs{i};
    [path, name, ext] = fileparts(file);
    
    switch method
        case 'bim2gipl'
            diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' name '_bim2gipl.log']));
            write_gipl(file, confs);
        case 'fix2gipl'
            diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' name '_fix2gipl.log']));
            write_gipl(file, confs);            
        case 'gipl2bim'
            diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' name '_gipl2bim.log']));            
            [bim, origin, vxsize] = read_gipl(file);
            origin = origin(1:3);
            new_name = [confs.OutDirectory '/' name '_bim.mat'];
            if exist(new_name,'file')
                prompt = {'Enter new filename:'};
                dlg_title = 'New File Name';
                num_lines = 1;
                def = {new_name};
                answer = inputdlg(prompt,dlg_title,num_lines,def);    
                new_name = answer{1};
            end
            save(new_name, 'bim', 'origin', 'vxsize');
        case 'smo2surf_para_meta'
            diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' name '_smo2surf.log']));
            write_surfNpara(file, confs);
        case 'surf_para_meta2smo'
            diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' name '_surf2smo.log']));

            if strcmp(name(end-3:end),'surf')
                para=sprintf('%s/%spara.meta',path,name(1:end-4));
                if ~exist(file,'file') | ~exist(para,'file')
                    disp('_para or _surf file does not exist');
                else
                    [vertices,faces] = read_meta(file);
                    [sph_verts, sph_faces] = read_meta(para);
                    new_name = [confs.OutDirectory '/' name '_smo.mat'];
                    if exist(new_name,'file')
                        prompt = {'Enter new filename:'};
                        dlg_title = 'New File Name';
                        num_lines = 1;
                        def = {new_name};
                        answer = inputdlg(prompt,dlg_title,num_lines,def);    
                        new_name = answer{1};
                    end
                    save(new_name,'vertices','sph_verts','faces');
                end
            end
        case 'des2meta_coef'
            answer = questdlg('SPHARM coefficients of negative order will be lost. Continue?', ...
                'Format Conversion: des2meta_coef','Proceed','Cancel','Proceed');
            diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' name '_des2meta.log']));            
            if strcmpi(answer,'Proceed')
                write_metaNcoef(file, confs);
            end
        case 'meta_coef2des'
            diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' name '_meta2des.log']));
            if strcmp(ext,'.meta')
                coef=sprintf('%s/%s.coef',path,name);
                if ~exist(file,'file') | ~exist(coef,'file')
                    disp('.meta or .coef file does not exist');
                else
                    [vertices,faces] = read_meta(file);
                    [fvec] = read_coef(coef);
                    new_name = [confs.OutDirectory '/' name '_PDM_des.mat'];
                    if exist(new_name,'file')
                        prompt = {'Enter new filename:'};
                        dlg_title = 'New File Name';
                        num_lines = 1;
                        def = {new_name};
                        answer = inputdlg(prompt,dlg_title,num_lines,def);    
                        new_name = answer{1};
                    end
                    sph_verts=[];
                    save(new_name,'vertices','fvec','faces','sph_verts');
                end
            end
        case 'ellalign_meta_coef2reg'
            diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' name '_meta2reg.log']));
            if strcmp(ext,'.meta')
                coef=sprintf('%s/%s.coef',path,name);
                if ~exist(file,'file') | ~exist(coef,'file')
                    disp('.meta or .coef file does not exist');
                else
                    [vertices,faces] = read_meta(file);
                    [fvec] = read_coef(coef);
                    new_name = [confs.OutDirectory '/' name '_PDM_reg.mat'];
                    if exist(new_name,'file')
                        prompt = {'Enter new filename:'};
                        dlg_title = 'New File Name';
                        num_lines = 1;
                        def = {new_name};
                        answer = inputdlg(prompt,dlg_title,num_lines,def);    
                        new_name = answer{1};
                    end
                    sph_verts=[];
                    save(new_name,'vertices','fvec','faces','sph_verts');
                end
            end
        case 'reg2procalign_meta'
            diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' name '_reg2meta.log']));
            write_procalign(file, confs);
    end
    diary('off');
    waitbar(i/numSbj)
end
close(h);

return;
