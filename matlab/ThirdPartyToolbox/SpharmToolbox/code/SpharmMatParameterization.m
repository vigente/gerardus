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

function outNames = SpharmMatParameterization(confs, objs, method)

numSbj = length(objs);
if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end
if (~exist([confs.OutDirectory '/Logs'],'dir'))
    mkdir([confs.OutDirectory '/Logs']);
end

outNames = {};

h = waitbar(0,'Please wait...');
for i = 1:numSbj
    file = objs{i};
    [path, name, ext] = fileparts(file);
    
    switch method
        case 'ParamCALD'
            diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' name '_ParamCALD.log']));
            % Call CALD method (Intial parameterization and smoothing)
            [vertices_new, faces_new, sph_verts_new, outNames{end+1}] = SpharmMatParamTriaMesh(file, confs);
            clear('vertices_new', 'faces_new', 'sph_verts_new');
        case 'ParamPDM'
            % Call SPHARM-PDM method
            diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' name '_ParamPDM.log']));
            optionSTR = '';
            for j = 1:length(confs.vars)
                vals = eval(sprintf('confs.%s', confs.vars{j}));
                if ~isempty(vals) & (confs.args(j) < 10)
                    vals = num2str(vals);
                    optionSTR = [optionSTR ' ' sprintf('-%s ', confs.vars{j}) vals];
                elseif (confs.args(j) < 200) & (confs.args(j) >= 10) & ~isempty(deblank(vals)) & ~strcmp(confs.vars{j}, 'others')
                    optionSTR = [optionSTR ' ' sprintf('-%s ', confs.vars{j}) vals];                    
                elseif (confs.args(j) < 200) & (confs.args(j) >= 10) & ~isempty(deblank(vals)) & strcmp(confs.vars{j}, 'others')
                    optionSTR = [optionSTR ' ' vals];                    
                    
                end
            end
            if ispc
                cSTR = sprintf('!%s %s %s;', confs.command, file, optionSTR);
                wdir = sprintf('%s/tmp',confs.path);
                if ~exist(wdir,'dir')
                    mkdir(wdir);
                end
                copyfile(file,wdir);
                [p,fname,ext,ver] = fileparts(file);
                cSTR = sprintf('!%s tmp/%s %s', confs.command, [fname ext], optionSTR);               
                new_name = [ 'tmp/' name '*.meta'];                
            elseif isunix
                new_name = [path '/' name '*.meta'];
                cSTR = sprintf('!./%s %s %s;', confs.command, file, optionSTR);
            elseif ismac
                new_name = [path '/' name '*.meta'];
                cSTR = sprintf('!./%s %s %s;', confs.command, file, optionSTR);
            end
            
            curPath = pwd;
            cd(confs.path);
            eval(cSTR);
            cd(curPath);
            if ispc
                movefile([confs.path '/' new_name],confs.OutDirectory);
            elseif isunix
                movefile(new_name,confs.OutDirectory);                
            end
            outNames{end+1} = [confs.OutDirectory '/' name '_surf.meta'];
    end
    diary('off');
    waitbar(i/numSbj)
end
close(h);

return