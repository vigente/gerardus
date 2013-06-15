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

disp('Ex0303: Perform Voxel Surfaces (PDM)')

SpharmMatDir = '../';
if ~ispc
    curP = pwd;
    cd(SpharmMatDir);
    SpharmMatDir = [pwd '/'];
    cd(curP);
end

CodeDir = [SpharmMatDir 'code'];
if exist([CodeDir '/sysConfig.mat'])
    load([CodeDir '/sysConfig.mat']);
else
    disp('PDM_path variable (sysConfig.mat) does not exist in SPHARM-MAT directory.')
    disp('Please set the path to SHARM-PDM for SPHARM-MAT and use this script.');
    return;
end
addpath(CodeDir);

DataFolder = 'data';

%% (1) Input data directory
InDataDir = [SpharmMatDir DataFolder '/Ex0303/hip03_gipl_fix'];


%% (2) List input file names
inFiles = dir([InDataDir '/*_fix.gipl']); inNames={};
for i=1:length(inFiles)
    inNames{end+1} = [InDataDir '/' inFiles(i).name];
end


%% (3) Perform PDM Parameterization for Voxel Surfaces
ParamPDM.command = 'GenParaMesh';
if ispc
    PDM_path = PDM_path_pc;
elseif isunix
    PDM_path = PDM_path_unix;
end

com = [PDM_path '/' ParamPDM.command]; 
ParamPDM.exist = 0;
if ispc
    if exist([com '.exe'],'file')
        ParamPDM.exist = 1;
    end
else
    if exist(com, 'file')
        ParamPDM.exist = 1;
    end
end

confs.vars = {'iter', 'label', 'OutDirectory', 'others'};
confs.args = [1 1 200 20];
confs.iter = 500;
confs.label = 1;
confs.OutDirectory = [SpharmMatDir DataFolder '/Ex0303/hip04_meta'];
confs.others = '';
confs.command = ParamPDM.command;
confs.path = PDM_path;

if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end

if (ParamPDM.exist == 0)
    errordlg(sprintf('%s function does not exist in the PDM directory, %s',PDM_Fix.command,PDM_path));
    rmpath(Codedir);
    return;
else
    SpharmMatParameterization(confs, inNames, 'ParamPDM');
end


%% (4) Format Conversion surf_para_meta2smo
inFiles = dir([confs.OutDirectory '/*.meta']); inNames={};
for i=1:length(inFiles)
    inNames{end+1} = [confs.OutDirectory '/' inFiles(i).name];
end

confs.OutDirectory = [SpharmMatDir DataFolder '/Ex0303/hip05_smo'];

if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end

SpharmMatUtilFormatConvert(confs, inNames, 'surf_para_meta2smo');


%% (5) Display output objects (optional)
inFiles = dir([confs.OutDirectory '/*_smo.mat']); inNames={};
for i=1:length(inFiles)
    inNames{end+1} = [confs.OutDirectory '/' inFiles(i).name];
end

    %Available values for Space- 'object';'param';'both'
dispConfs.Space = 'both';
    % Available values for Mesh- 'orig';'quad32';'quad64';'quad128';'quad256';'quad512';'icosa1';'icosa2'; ...
    %    'icosa3';'icosa4';'icosa5';'icosa6'
dispConfs.Mesh = 'orig';  
    % Available values for Shape- 'solid';'mesh';'both'
dispConfs.Shade = 'both';
    % Available values for Overlay- 'none';'adc_paramap'
dispConfs.Overlay = 'adc_paramap';
    % Available values for Export- 'screen';'png';'both'
dispConfs.Export = 'png';
dispConfs.Degree = [];
dispConfs.Template = '';

SpharmMatUtilDisplayObjs(dispConfs, inNames, CodeDir);


rmpath(CodeDir);

clear;