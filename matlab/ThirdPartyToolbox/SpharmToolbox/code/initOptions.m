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

function userdataHandle = initOptions(userdataHandle)

% Start initialization
% Not implemented: 'ParamQuad', 'ParamConf'
userdataHandle.keyword.list = {'ParamCALD', 'ParamPDM', 'ParamQuad', 'ParamConf', 'ExpLSF', 'ExpIRF', ...
    'AligSHREC', 'AligFOE', 'AligLandmark', 'ExpAlig', 't_map', 'PCA', 't_map+PCA', ... 
    'PDM_stat', 'Import', 'DisplayRes', 'bim2gipl', 'fix2gipl', 'gipl2bim', 'smo2surf_para_meta', ...
    'surf_para_meta2smo', 'des2meta_coef', 'meta_coef2des', 'reg2procalign_meta', 'InHouse_Fix', ...
    'PDM_Fix', 'Config'};
userdataHandle.keyword.selected = '';

userdataHandle.mesh={'orig';'quad32';'quad64';'quad128';'quad256';'quad512';'icosa1';'icosa2'; ...
    'icosa3';'icosa4';'icosa5';'icosa6'};
userdataHandle.space={'object';'param';'both'};
userdataHandle.shade={'solid';'mesh';'both'};
userdataHandle.overlay={'none';'adc_paramap'};
userdataHandle.colormap={'jet';'hot';'summer';'cool';'autumn';'winter';'spring'};
userdataHandle.export={'screen';'png';'both'};

% CALD Parameterization
userdataHandle.ParamCALD.vars = {'MeshGridSize','MaxSPHARMDegree','Tolerance','Smoothing','Iteration', ... 
        'LocalIteration','t_major','SelectDiagonal','OutDirectory'};
userdataHandle.ParamCALD.args = [1 1 1 1 1 1 10 10 200];  % Encode the number and type of each variables  
userdataHandle.ParamCALD.inFilter = {...
    '*_bim.mat;*_fix.mat;*_obj.mat', 'Binary and Mesh Objects (*_bim.mat;*_fix.mat;*_obj.mat)'; ...
    '*_obj.mat','Mesh Objects (*_obj.mat)'; ...
    '*_bim.mat;*_fix.mat', 'Bianry Objects (*_bim.mat, *_fix.mat)'; ...
    '*.mat', 'All Objects (*.mat)' ...
    };
userdataHandle.ParamCALD.default={'50', '6', '2', '2', '100', '10', '', '','./'};

for i=1:length(userdataHandle.ParamCALD.vars)
    varName = sprintf('userdataHandle.ParamCALD.%s', userdataHandle.ParamCALD.vars{i});
    if userdataHandle.ParamCALD.args(i) < 10
        strCC = sprintf('%s=str2num(userdataHandle.ParamCALD.default{%d});', varName, i);
        eval(strCC);
    else
        strCC = sprintf('%s=(userdataHandle.ParamCALD.default{%d});', varName, i);
        eval(strCC);
    end
end

% PDM Parameterization (GenParamMesh)
userdataHandle.ParamPDM.vars = {'iter', 'label', 'OutDirectory', 'others'};
userdataHandle.ParamPDM.args = [1 1 200 20];
userdataHandle.ParamPDM.inFilter = {'*.gipl;*.hdr', 'Bianry Objects (*.gipl, *hdr)'};
userdataHandle.ParamPDM.default = {'500','1','./',''};
userdataHandle.ParamPDM.path = '';
userdataHandle.ParamPDM.command = 'GenParaMesh';
userdataHandle.ParamPDM.exist = 0;

for i=1:length(userdataHandle.ParamPDM.vars)
    varName = sprintf('userdataHandle.ParamPDM.%s', userdataHandle.ParamPDM.vars{i});
    if userdataHandle.ParamPDM.args(i) < 10
        strCC = sprintf('%s=str2num(userdataHandle.ParamPDM.default{%d});', varName, i);
        eval(strCC);
    else
        strCC = sprintf('%s=(userdataHandle.ParamPDM.default{%d});', varName, i);
        eval(strCC);
    end
end

% Parameterization using a similar method to Brechbuhler's methd
userdataHandle.ParamQuad.vars = {'NRS_step','OutDirectory'};
userdataHandle.ParamQuad.args = [1 200];  % Encode the number and type of each variables  
userdataHandle.ParamQuad.inFilter = {'*_bim.mat;*_fix.mat', 'Bianry Objects (*_bim.mat, *_fix.mat)'; ...
   '*_obj.mat','Mesh Objects (*_obj.mat)'; '*.mat', 'All Objects (*.mat)'};
userdataHandle.ParamQuad.default={'2','./'};

for i=1:length(userdataHandle.ParamQuad.vars)
    varName = sprintf('userdataHandle.ParamQuad.%s', userdataHandle.ParamQuad.vars{i});
    if userdataHandle.ParamQuad.args(i) < 10
        strCC = sprintf('%s=str2num(userdataHandle.ParamQuad.default{%d});', varName, i);
        eval(strCC);
    else
        strCC = sprintf('%s=(userdataHandle.ParamQuad.default{%d});', varName, i);
        eval(strCC);
    end
end

% Expansion in the sense of least mean square
userdataHandle.ExpLSF.vars = {'MaxSPHARMDegree','OutDirectory'};
userdataHandle.ExpLSF.args = [1 200];
userdataHandle.ExpLSF.inFilter = {'*_ini.mat;*_smo.mat', 'Parameterized Objects (*_ini.mat, *_smo.mat)'; '*.mat', 'All Objects (*.mat)'};
userdataHandle.ExpLSF.default = {'15', './'};

for i=1:length(userdataHandle.ExpLSF.vars)
    varName = sprintf('userdataHandle.ExpLSF.%s', userdataHandle.ExpLSF.vars{i});
    if userdataHandle.ExpLSF.args(i) < 10
        strCC = sprintf('%s=str2num(userdataHandle.ExpLSF.default{%d});', varName, i);
        eval(strCC);
    else
        strCC = sprintf('%s=(userdataHandle.ExpLSF.default{%d});', varName, i);
        eval(strCC);
    end
end

% SHREC Alignment
userdataHandle.AligSHREC.vars = {'Template', 'MaxSPHARMDegree', 'GroupAlpha', 'NormalizeSize', ...
    'BaseRes', 'HierarchyStep', 'HierarchyDepth', 'Top_K', 'GammaRes', 'OutDirectory'};
userdataHandle.AligSHREC.args = [100 1 1 10 1 1 1 1 1 200];
userdataHandle.AligSHREC.inFilter = {'*_des.mat;*_prm.mat', 'SPHARM Descriptors (*_des.mat,*_prm.mat)'; '*.mat', 'All Objects (*.mat)'};
userdataHandle.AligSHREC.default = {'','12','100','','1','1','3','1','2','./'};

for i=1:length(userdataHandle.AligSHREC.vars)
    varName = sprintf('userdataHandle.AligSHREC.%s', userdataHandle.AligSHREC.vars{i});
    if userdataHandle.AligSHREC.args(i) < 10
        strCC = sprintf('%s=str2num(userdataHandle.AligSHREC.default{%d});', varName, i);
        eval(strCC);
    else
        strCC = sprintf('%s=(userdataHandle.AligSHREC.default{%d});', varName, i);
        eval(strCC);
    end
end

% First order ellipsoid alignment
userdataHandle.AligFOE.vars = {'CPoint', 'NPole', 'MaxSPHARMDegree', 'OutDirectory'};
userdataHandle.AligFOE.args = [10 10 1 200];
userdataHandle.AligFOE.inFilter = {'*_des.mat', 'SPHARM Descriptors (*_des.mat)'; '*.mat', 'All Objects (*.mat)'};
userdataHandle.AligFOE.default = {'', '', '12', './'};

for i=1:length(userdataHandle.AligFOE.vars)
    varName = sprintf('userdataHandle.AligFOE.%s', userdataHandle.AligFOE.vars{i});
    if userdataHandle.AligFOE.args(i) < 10
        strCC = sprintf('%s=str2num(userdataHandle.AligFOE.default{%d});', varName, i);
        eval(strCC);
    else
        strCC = sprintf('%s=(userdataHandle.AligFOE.default{%d});', varName, i);
        eval(strCC);
    end
end

% PDM expansion and alignment (ParaToSPHARMMesh)
userdataHandle.ExpAlig.vars = {'flipTemplate', 'subdivLevel', 'spharmDegree', ...
        'regTemplate', 'FinalFlip', 'OutDirectory', 'others'};
userdataHandle.ExpAlig.args = [100 1 1 100 1 200 20];
userdataHandle.ExpAlig.inFilter = {'*_surf.meta', 'Parameterized Objects (*_surf.meta)'};
userdataHandle.ExpAlig.default = {'','20','12','','0','./',''};    
userdataHandle.ExpAlig.path = '';
userdataHandle.ExpAlig.command = 'ParaToSPHARMMesh';
userdataHandle.ExpAlig.exist = 0;

for i=1:length(userdataHandle.ExpAlig.vars)
    varName = sprintf('userdataHandle.ExpAlig.%s', userdataHandle.ExpAlig.vars{i});
    if userdataHandle.ExpAlig.args(i) < 10
        strCC = sprintf('%s=str2num(userdataHandle.ExpAlig.default{%d});', varName, i);
        eval(strCC);
    else
        strCC = sprintf('%s=(userdataHandle.ExpAlig.default{%d});', varName, i);
        eval(strCC);
    end
end

% t-statistical analysis 
userdataHandle.t_map.vars = {'Atlas', 'Smoothing_FWHM', 'EqualVariance', 'Signal', 'SampleMesh', 'OutputNamePrefix','OutDirectory','GroupIDs'};
userdataHandle.t_map.args = [100 1 10 10 10 10 200 10];
userdataHandle.t_map.inFilter = {'*.mat', 'Registered Objects (*.mat)'; '*.mat', 'All Objects (*.mat)'};
userdataHandle.t_map.default = {'./atlas.mat','5','','','','t_map','./','Ctrl,PT'};

for i=1:length(userdataHandle.t_map.vars)
    varName = sprintf('userdataHandle.t_map.%s', userdataHandle.t_map.vars{i});
    if userdataHandle.t_map.args(i) < 10
        strCC = sprintf('%s=str2num(userdataHandle.t_map.default{%d});', varName, i);
        eval(strCC);
    else
        strCC = sprintf('%s=(userdataHandle.t_map.default{%d});', varName, i);
        eval(strCC);
    end
end

% PCA analysis 
userdataHandle.PCA.vars = {'GroupID','OutputName', 'OutDirectory'};
userdataHandle.PCA.args = [10 10 200];
userdataHandle.PCA.inFilter = {'*.mat', 'Registered Objects (*.mat)'; '*.mat', 'All Objects (*.mat)'};
userdataHandle.PCA.default = {'All','PCA_stat.mat','./'};

for i=1:length(userdataHandle.PCA.vars)
    varName = sprintf('userdataHandle.PCA.%s', userdataHandle.PCA.vars{i});
    if userdataHandle.PCA.args(i) < 10
        strCC = sprintf('%s=str2num(userdataHandle.PCA.default{%d});', varName, i);
        eval(strCC);
    else
        strCC = sprintf('%s=(userdataHandle.PCA.default{%d});', varName, i);
        eval(strCC);
    end
end

% Display T-stat maps
userdataHandle.res_t_map.vars = {'Threshold_p_value', 'Overlay','Colormap'};
userdataHandle.res_t_map.args = [1 10 10];
userdataHandle.res_t_map.inFilter = {'t_map_*.mat',  'T-Stat Results (t_map_*.mat)'; ... 
    '*.mat', 'All Objects (*.mat)'};
userdataHandle.res_t_map.default = {'0.05', '', ''};

for i=1:length(userdataHandle.res_t_map.vars)
    varName = sprintf('userdataHandle.res_t_map.%s', userdataHandle.res_t_map.vars{i});
    if userdataHandle.res_t_map.args(i) < 10
        strCC = sprintf('%s=str2num(userdataHandle.res_t_map.default{%d});', varName, i);
        eval(strCC);
    else
        strCC = sprintf('%s=(userdataHandle.res_t_map.default{%d});', varName, i);
        eval(strCC);
    end
end

% Display PCA results
userdataHandle.res_PCA.vars = {'Level','Sigma','Mesh','MaxSPHARMDegree'};
userdataHandle.res_PCA.args = [1 1 10 1];
userdataHandle.res_PCA.inFilter = {'PCA_*.mat',  'PCA Results (PCA_*.mat)'; ... 
    '*.mat', 'All Objects (*.mat)'};
userdataHandle.res_PCA.default = {'4','3','','15'};

for i=1:length(userdataHandle.res_PCA.vars)
    varName = sprintf('userdataHandle.res_PCA.%s', userdataHandle.res_PCA.vars{i});
    if userdataHandle.res_PCA.args(i) < 10
        strCC = sprintf('%s=str2num(userdataHandle.res_PCA.default{%d});', varName, i);
        eval(strCC);
    else
        strCC = sprintf('%s=(userdataHandle.res_PCA.default{%d});', varName, i);
        eval(strCC);
    end
end

% PDM topology fix (SegPostProcess)
userdataHandle.PDM_Fix.vars = {'space','OutDirectory', 'others'};
userdataHandle.PDM_Fix.args = [3 200 20];
userdataHandle.PDM_Fix.inFilter = {'*.gipl;*.hdr', 'Bianry Objects (*.gipl, *hdr)'};
userdataHandle.PDM_Fix.default = {'0.75,0.75,0.75','./',''};
userdataHandle.PDM_Fix.path = '';
userdataHandle.PDM_Fix.command = 'SegPostProcess';
userdataHandle.PDM_Fix.exist = 0;

for i=1:length(userdataHandle.PDM_Fix.vars)
    varName = sprintf('userdataHandle.PDM_Fix.%s', userdataHandle.PDM_Fix.vars{i});
    if userdataHandle.PDM_Fix.args(i) < 10
        strCC = sprintf('%s=str2num(userdataHandle.PDM_Fix.default{%d});', varName, i);
        eval(strCC);
    else
        strCC = sprintf('%s=(userdataHandle.PDM_Fix.default{%d});', varName, i);
        eval(strCC);
    end
end

% InHouse topology fixing method
userdataHandle.InHouse_Fix.vars = {'Connectivity','Epsilon','OutDirectory'};
userdataHandle.InHouse_Fix.args = [10 1 200];
userdataHandle.InHouse_Fix.inFilter = {'*_bim.mat;*_fix.mat', 'Bianry Objects (*_bim.mat, *_fix.mat)'; ...
    '*.mat', 'All Objects (*.mat)'};
userdataHandle.InHouse_Fix.default = {'','1.5','./'};

for i=1:length(userdataHandle.InHouse_Fix.vars)
    varName = sprintf('userdataHandle.InHouse_Fix.%s', userdataHandle.InHouse_Fix.vars{i});
    if userdataHandle.InHouse_Fix.args(i) < 10
        strCC = sprintf('%s=str2num(userdataHandle.InHouse_Fix.default{%d});', varName, i);
        eval(strCC);
    else
        strCC = sprintf('%s=(userdataHandle.InHouse_Fix.default{%d});', varName, i);
        eval(strCC);
    end
end

% System configuration    
userdataHandle.Config.vars = {'PDM_path_pc','PDM_path_unix'};
userdataHandle.Config.args = [200,200];
userdataHandle.Config.default = {'./','./'};
userdataHandle.Config.filename = 'sysConfig.mat';

for i=1:length(userdataHandle.Config.vars)
    varName = sprintf('userdataHandle.Config.%s', userdataHandle.Config.vars{i});
    if userdataHandle.Config.args(i) < 10
        strCC = sprintf('%s=str2num(userdataHandle.Config.default{%d});', varName, i);
        eval(strCC);
    else
        strCC = sprintf('%s=(userdataHandle.Config.default{%d});', varName, i);
        eval(strCC);
    end
end

userdataHandle.Config.loaded = 0;

if exist(fullfile(pwd, userdataHandle.Config.filename), 'file')
    load(fullfile(pwd, userdataHandle.Config.filename));
    if ispc
        PDM_path = PDM_path_pc;
    elseif isunix
        PDM_path = PDM_path_unix;        
    end

    userdataHandle.Config.PDM_path_pc = PDM_path_pc;    
    userdataHandle.Config.PDM_path_unix = PDM_path_unix;    
    
    if exist('PDM_path','var')
        userdataHandle.Config.loaded = 1;
        userdataHandle.Config.PDM_path = PDM_path;
        userdataHandle.ParamPDM.path = PDM_path;
        userdataHandle.PDM_Fix.path = PDM_path;
        userdataHandle.ExpAlig.path = PDM_path;
%         userdataHandle.PDM_stat.path = PDM_path;
    end        
    
    com1 = [PDM_path '/' userdataHandle.ParamPDM.command];
    com2 = [PDM_path '/' userdataHandle.ExpAlig.command];    
%     com3 = [PDM_path '/' userdataHandle.PDM_stat.command];
    com4 = [PDM_path '/' userdataHandle.PDM_Fix.command]; 
    
    if ispc
        com1 = [PDM_path '/' userdataHandle.ParamPDM.command '.exe'];
        com2 = [PDM_path '/' userdataHandle.ExpAlig.command '.exe'];    
        com4 = [PDM_path '/' userdataHandle.PDM_Fix.command '.exe']; 
    end
    if exist(com1, 'file')
        userdataHandle.ParamPDM.exist =1;
    else
        userdataHandle.ParamPDM.exist =0;
    end
    if exist(com2, 'file')
        userdataHandle.ExpAlig.exist =1;
    else
        userdataHandle.ExpAlig.exist =0;
    end
    if exist(com4, 'file')
        userdataHandle.PDM_Fix.exist =1;
    else
        userdataHandle.PDM_Fix.exist =0;
    end
end

% Display objects
userdataHandle.DisplayObjs.vars = {'Space','Mesh','Shade','Overlay','Export','Degree','Template'};
userdataHandle.DisplayObjs.args = [10 10 10 10 10 1 100];
userdataHandle.DisplayObjs.inFilter = {...
    '*.mat', 'All Objects (*.mat)'; ...
    '*_bim.mat; *_fix.mat', 'Bianry Objects (*_bim.mat, *_fix.mat)'; ...
    '*_obj.mat', 'Mesh Objects (*_obj.mat)'; ...
    '*_ini.mat; *_smo.mat', 'Parameterized Objects (*_ini.mat, *_smo.mat)'; ...
    '*_des.mat', 'SPHARM Descriptors (*_des.mat)'; ...
    '*_reg.mat; *_nor.mat',  'Registered Objects (*_reg.mat, *_nor.mat)' ... 
    };
userdataHandle.DisplayObjs.default = {'','','','','','',''};

for i=1:length(userdataHandle.DisplayObjs.vars)
    varName = sprintf('userdataHandle.DisplayObjs.%s', userdataHandle.DisplayObjs.vars{i});
    if userdataHandle.DisplayObjs.args(i) < 10
        strCC = sprintf('%s=str2num(userdataHandle.DisplayObjs.default{%d});', varName, i);
        eval(strCC);
    else
        strCC = sprintf('%s=(userdataHandle.DisplayObjs.default{%d});', varName, i);
        eval(strCC);
    end
end

%  Scale objects
userdataHandle.ScaleObjs.vars = {'ScalingFactor','OutDirectory'};
userdataHandle.ScaleObjs.args = [100 200];
userdataHandle.ScaleObjs.inFilter = {'*_des.mat; *_prm.mat; *_reg.mat', 'SPHARM Descriptors (*_des.mat, *_prm.mat, *_reg.mat)'; ...
    '*_obj.mat; *_ini.mat; *_smo.mat',  'Meshed Objects (*_obj.mat, *_ini.mat, *_smo.mat)'; ... 
    '*.mat', 'All Objects (*.mat)'};
userdataHandle.ScaleObjs.default = {'','./'};

for i=1:length(userdataHandle.ScaleObjs.vars)
    varName = sprintf('userdataHandle.ScaleObjs.%s', userdataHandle.ScaleObjs.vars{i});
    if userdataHandle.ScaleObjs.args(i) < 10
        strCC = sprintf('%s=str2num(userdataHandle.ScaleObjs.default{%d});', varName, i);
        eval(strCC);
    else
        strCC = sprintf('%s=(userdataHandle.ScaleObjs.default{%d});', varName, i);
        eval(strCC);
    end
end

% Create average object
userdataHandle.AverageObjs.vars = {'OutputName','OutDirectory'};
userdataHandle.AverageObjs.args = [10 200];
userdataHandle.AverageObjs.inFilter = {'*_des.mat; *_prm.mat; *_reg.mat', 'SPHARM Descriptors (*_des.mat, *_prm.mat, *_reg.mat)'; ...
    '*.mat', 'All Objects (*.mat)'};
userdataHandle.AverageObjs.default = {'atlas.mat','./'};

for i=1:length(userdataHandle.AverageObjs.vars)
    varName = sprintf('userdataHandle.AverageObjs.%s', userdataHandle.AverageObjs.vars{i});
    if userdataHandle.AverageObjs.args(i) < 10
        strCC = sprintf('%s=str2num(userdataHandle.AverageObjs.default{%d});', varName, i);
        eval(strCC);
    else
        strCC = sprintf('%s=(userdataHandle.AverageObjs.default{%d});', varName, i);
        eval(strCC);
    end
end

% Import objects
userdataHandle.Import.vars = {'ResampleFactor','OutDirectory'};
userdataHandle.Import.args = [1 200];
userdataHandle.Import.inFilter = {'*.stl', 'Mesh Objects (*.stl)'; ...
    '*.m', 'Surface-mesh Objects (*.m)'; ...
    '*.hdr; *.nii', 'Binary Analyze/NIfTI Objects (*.hdr, *.nii)'};% ...
%    '*.stl; *.m; *.hdr; *.nii', 'All Eligible Files (*.stl, *.m, *.hdr, *.nii)'};
userdataHandle.Import.default = {'1','./'};

for i=1:length(userdataHandle.Import.vars)
    varName = sprintf('userdataHandle.Import.%s', userdataHandle.Import.vars{i});
    if userdataHandle.Import.args(i) < 10
        strCC = sprintf('%s=str2num(userdataHandle.Import.default{%d});', varName, i);
        eval(strCC);
    else
        strCC = sprintf('%s=(userdataHandle.Import.default{%d});', varName, i);
        eval(strCC);
    end
end

% Convert file formats
userdataHandle.bim2gipl.inFilter = {'*_bim.mat', 'Bianry Objects (*_bim.mat)'; ...
    '*.mat', 'All Objects (*.mat)'};
userdataHandle.bim2gipl.args=[200];
userdataHandle.bim2gipl.vars = {'OutDirectory'};
userdataHandle.bim2gipl.OutDirectory = './';

userdataHandle.fix2gipl.inFilter = {'*_fix.mat', 'Bianry Objects (*_fix.mat)'; ...
    '*.mat', 'All Objects (*.mat)'};
userdataHandle.fix2gipl.args=[200];
userdataHandle.fix2gipl.vars = {'OutDirectory'};
userdataHandle.fix2gipl.OutDirectory = './';

userdataHandle.gipl2bim.inFilter = {'*.gipl', 'Bianry Objects (*.gipl)'};
userdataHandle.gipl2bim.args=[200];
userdataHandle.gipl2bim.vars = {'OutDirectory'};
userdataHandle.gipl2bim.OutDirectory = './';

userdataHandle.smo2surf_para_meta.inFilter = {'*_smo.mat', 'Parameterized Objects (*_smo.mat)'};
userdataHandle.smo2surf_para_meta.args=[200];
userdataHandle.smo2surf_para_meta.vars = {'OutDirectory'};
userdataHandle.smo2surf_para_meta.OutDirectory = './';

userdataHandle.surf_para_meta2smo.inFilter = {'*_surf.meta;*_para.meta', 'Parameterized Objects (*_surf.meta, *_para.meta)'};
userdataHandle.surf_para_meta2smo.args=[200];
userdataHandle.surf_para_meta2smo.vars = {'OutDirectory'};
userdataHandle.surf_para_meta2smo.OutDirectory = './';

userdataHandle.des2meta_coef.inFilter = {'*_des.mat', 'SPHARM Descriptors (*_des.mat)'};
userdataHandle.des2meta_coef.args=[200];
userdataHandle.des2meta_coef.vars = {'OutDirectory'};
userdataHandle.des2meta_coef.OutDirectory = './';

userdataHandle.meta_coef2des.inFilter = {'*_surfSPHARM.meta;*_surfSPHARM.coef', 'SPHARM Descriptors (*_surfSPHARM.meta, *_surfSPHARM.coef)'};
userdataHandle.meta_coef2des.args=[200];
userdataHandle.meta_coef2des.vars = {'OutDirectory'};
userdataHandle.meta_coef2des.OutDirectory = './';

userdataHandle.ellalign_meta_coef2reg.inFilter = {'*_ellalign.meta;*_ellalign.coef', 'Registered SPHARM Descriptors (*_ellalign.meta, *_ellalign.coef)'};
userdataHandle.ellalign_meta_coef2reg.args=[200];
userdataHandle.ellalign_meta_coef2reg.vars = {'OutDirectory'};
userdataHandle.ellalign_meta_coef2reg.OutDirectory = './';

% reg2procalign_meta option (mesh selection)
userdataHandle.reg2procalign_meta.vars = {'SampleMesh','OutDirectory'};
userdataHandle.reg2procalign_meta.args = [10,200];
userdataHandle.reg2procalign_meta.inFilter = {'*_reg.mat;*_nor.mat', 'Registered Objects (*_reg.mat, *_nor.mat)'};
userdataHandle.reg2procalign_meta.default = {'','./'};

for i=1:length(userdataHandle.reg2procalign_meta.vars)
    varName = sprintf('userdataHandle.reg2procalign_meta.%s', userdataHandle.reg2procalign_meta.vars{i});
    if userdataHandle.reg2procalign_meta.args(i) < 10
        strCC = sprintf('%s=str2num(userdataHandle.reg2procalign_meta.default{%d});', varName, i);
        eval(strCC);
    else
        strCC = sprintf('%s=(userdataHandle.reg2procalign_meta.default{%d});', varName, i);
        eval(strCC);
    end
end

userdataHandle.others.inFilter = {'*.txt', 'List of Objects (*.txt)'};  %'t_map' 'PCA' 't_map+PCA' 

return;