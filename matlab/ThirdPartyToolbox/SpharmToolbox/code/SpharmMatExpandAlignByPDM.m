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

function SpharmMatExpandAlignPDM(confs, objs)
numSbj = length(objs);
if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end
if (~exist([confs.OutDirectory '/Logs'],'dir'))
    mkdir([confs.OutDirectory '/Logs']);
end

h = waitbar(0,'Please wait...');
for i = 1:numSbj
    file1 = objs{i};
    [pa, name, ext, v] = fileparts(file1);
    file2 = [pa '/' name(1:end-4) 'para' ext];
    diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' name '_ExpPDM.log']));
    
    if ~exist(file2,'file')
        disp(sprintf('Corresponding parameter file to %s doesnot exist',name));
    else
    
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
            cSTR = sprintf('!%s %s %s %s', confs.command, file1, file2, optionSTR);
            disp('This function is not working on MS Windows. Please use Linux version.');
        elseif isunix
            preload = '/lib64/libgcc_s.so.1';
            if exist(preload, 'file')
                setenv('LD_PRELOAD', preload);
            else
                setenv('LD_PRELOAD', 'lib/libcc_s_so.1');
            end
            cSTR = sprintf('!./%s %s %s %s', confs.command, file1, file2, optionSTR);
            curPath = pwd;
            cd(confs.path);
                eval(cSTR);
            cd(curPath);
            movefile([pa '/' name 'SPHARM*.meta'],confs.OutDirectory);
            movefile([pa '/' name 'SPHARM*.coef'],confs.OutDirectory);  
        elseif ismac
            disp('Not implemented yet.');
        end
    end
    diary('off');
    
    waitbar(i/numSbj)
end
close(h);

return;
