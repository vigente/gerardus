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

function batchJobs(confs, objs, method)

numSbj = length(objs);
file1 = objs{1};
[path, name, ext, ver] = fileparts(file1);

if strcmp(method,'PatchPDM')
    if (strcmp('.gipl', ext) | strcmp('.hdr',ext)) & strcmp(confs.ParamPDM, 'None')
        disp('Binary input files can not be processed without PDM parameterization step.');
        disp('Please use individual processing steps in non-batching mode.');
        return;
    end

    if strcmp(ext,'.meta') & strcmp(confs.ExpPDM,'None')
        disp('Parameterized input files can not be processed without PDM expansion and alignment steps.');
        disp('Please use individual processing steps in non-batching mode.');        
        return;
        
    end
        
    if ~exist(confs.OutDirectory,'dir')
        mkdir(confs.OutDirectory);
    end
    if (~exist([confs.OutDirectory '/Logs'],'dir'))
        mkdir([confs.OutDirectory '/Logs']);
    end
else
    type = name(end-2:end);
    if (strcmp(type,'bim') | strcmp(type,'fix')) & strcmp(confs.ParamMethod, 'None')
        disp('Binary input files can not be processed without MAT parameterization step.');
        disp('Please use individual processing steps in non-batching mode.');        
        return;
    end

    if strcmp(type,'obj') & strcmp(confs.ParamMethod, 'None')
        disp('Meshed input files can not be processed without MAT parameterization step.');
        disp('Please use individual processing steps in non-batching mode.');        
        return;
    end
    
    if (strcmp(type,'ini') | strcmp(type,'smo')) & strcmp(confs.ExpMethod, 'None')
        disp('Smoothed input files can not be processed without MAT expansion step.');
        disp('Please use individual processing steps in non-batching mode.');        
        return;
    end
    
    if strcmp(type,'des') & strcmp(confs.SMOMethod, 'None')
        disp('SPHARM expanded input files can not be processed without MAT alignment step.');
        disp('Please use individual processing steps in non-batching mode.');        
        return;
    end

end

options.objs = objs;
options = initOptions(options);

switch method
    case 'BatchMAT'
        % Fix bad topology
        if (strcmp(type,'bim') | strcmp(type,'fix')) & strcmp(confs.TopologyFix, 'InHouse_Fix') & ~isempty(confs.TPOption)
            disp('Fixing bad topology');
            % Load options and run InHouse_Fix
            options.InHouse_Fix = loadOptions(confs.TPOption, options.InHouse_Fix);
            options.objs = topologyFix(options.InHouse_Fix, options.objs, 'InHouse_Fix');
            if ~isempty(options.objs), type = 'fix'; end;
        end

        % Parameterize objects
        if (strcmp(type,'bim') | strcmp(type,'fix') | strcmp(type,'obj')) & strcmp(confs.ParamMethod, 'ParamCALD') & ~isempty(confs.PMOption)
            disp('Parameterizing objects using CALD');
            % Load options and run ParamCALD
            options.ParamCALD = loadOptions(confs.PMOption, options.ParamCALD);
            options.objs = parameterizeObjs(options.ParamCALD, options.objs, 'ParamCALD');
            if ~isempty(options.objs), type = 'smo'; end;
        elseif (strcmp(type,'bim') | strcmp(type,'fix') | strcmp(type,'obj')) & strcmp(confs.ParamMethod, 'ParamQuad') & ~isempty(confs.PMOption)
            disp('Parameterizing objects using in-house quad');
            % Load options and run ParamQuad
            options.ParamQuad = loadOptions(confs.PMOption, options.ParamQuad);
            options.objs = parameterizeObjs(options.ParamQuad, options.objs, 'ParamQuad');
            if ~isempty(options.objs), type = 'smo'; end;
        end

        % SPHARM Expansion
        if (strcmp(type,'smo') | strcmp(type,'ini')) & strcmp(confs.ExpMethod, 'ExpLSF') & ~isempty(confs.ExpOption)
            disp('Calculating SPHARM coefficients');
            % Load options and run ExpLSF
            options.ExpLSF = loadOptions(confs.ExpOption, options.ExpLSF);
            options.objs = expandObjs(options.ExpLSF, options.objs, 'ExpLSF');
            if ~isempty(options.objs), type = 'des'; end;
        end

        % Alignment
        if strcmp(type,'des') & strcmp(confs.SMOMethod, 'AligSHREC') & ~isempty(confs.SMOOption)
            disp('Aligning objects using SHREC');
            % Load options and run AligSHREC
            options.AligSHREC = loadOptions(confs.SMOOption, options.AligSHREC);
            options.objs = alignObjs(options.AligSHREC, options.objs, 'AligSHREC');
            if ~isempty(options.objs), type = 'reg'; end;
        elseif strcmp(type,'des') & strcmp(confs.SMOMethod, 'AligFOE') & ~isempty(confs.SMOOption)
            disp('Aligning objects using FOE');
            % Load options and run AligFOE
            options.AligFOE = loadOptions(confs.SMOOption, options.AligFOE);
            options.objs = alignObjs(options.AligFOE, options.objs, 'AligFOE');
            if ~isempty(options.objs), type = 'reg'; end;
        end
        
    case 'BatchPDM'
        if (options.Config.loaded == 0)
            errordlg('Path to PDM program is not set.\n You need to set the path to PDM program in config menu of tool popup menu.');
            return;
        end

        % Fix bad topology
        if (options.PDM_Fix.exist == 0)
            errordlg(sprintf('%s function does not exist in the PDM directory, %s',options.PDM_Fix.command, options.PDM_path));
            return;
        elseif (strcmp('.gipl', ext) | strcmp('.hdr',ext)) & strcmp(confs.TopologyFix, 'SegPostProcess')
            disp('Fixing bad topology');
            % Load options and run PDM_Fix
            options.PDM_Fix.space=[];
            options.PDM_Fix.OutDirectory = confs.OutDirectory;
            options.PDM_Fix.others= confs.TPOption;
            options.objs = topologyFix(options.PDM_Fix, options.objs, 'PDM_Fix');
        end

        % Parameterize objects
        if (options.ParamPDM.exist == 0)
            errordlg(sprintf('%s function does not exist in the PDM directory, %s',options.ParamPDM.command, options.PDM_path));
            return;
        elseif isempty(options.objs)
            errordlg('Error during topology fixing. List of filenames is empty');
            return;
        elseif (strcmp('.gipl', ext) | strcmp('.hdr',ext)) & strcmp(confs.ParamPDM, 'GenParaMesh')
            disp('Parameterize objects');
            % Load options and run ParamPDM
            options.ParamPDM.iter=[];
            options.ParamPDM.label = [];
            options.ParamPDM.OutDirectory = confs.OutDirectory;
            options.ParamPDM.others= confs.ParamOption;
            options.objs = parameterizeObjs(options.ParamPDM, options.objs, 'ParamPDM');
            if ~isempty(optionis.objs), [pa,na,ext,ver] = fileparts(options.objs{1}); end;
        end

        % SPHARM Expansion and alignment
        if (options.ExpAlig.exist == 0)
            errordlg(sprintf('%s function does not exist in the PDM directory, %s',options.ExpAlig.command, options.Config.PDM_path));
        elseif isempty(options.objs)
            errordlg('Error during parameterization. List of filenames is empty');
            return;
        elseif strcmp('_surf.meta', [na(end-4:end) ext]) & strcmp(confs.ExpPDM, 'ParaToSPHARMMesh')
            disp('Calculate SPHARM coefficients and align objects');
            options.flipTemplate = '';
            options.subdivLevel = [];
            options.spharmDegree = [];
            options.regTemplate = '';
            options.FinalFlip = [];
            options.ExpAlig.OutDirectory = confs.OutDirectory;
            options.ExpAlig.others = confs.ExpOption;
    
            expandObjs_PDM(options.ExpAlig, options.objs);
        end
end

clear('options');

return;
