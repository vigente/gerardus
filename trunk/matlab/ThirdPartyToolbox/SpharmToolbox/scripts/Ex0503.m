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

disp('Ex0503: Perform SPHARM-PDM FOE alignment')

if ispc
    disp('SPHARM-PDM Alignment does not work on Windows');
    return;
end

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
InDataDir = [SpharmMatDir DataFolder '/Ex0503/hip04_meta'];


%% (2) List input file names
inFiles = dir([InDataDir '/*_surf.meta']); inNames={};
for i=1:length(inFiles)
    inNames{end+1} = [InDataDir '/' inFiles(i).name];
end


%% (3) Perform SPHARM-PDM FOE alignment
ExpAlig.command = 'ParaToSPHARMMesh';
if ispc
    PDM_path = PDM_path_pc;
elseif isunix
    PDM_path = PDM_path_unix;
end

com = [PDM_path '/' ExpAlig.command]; 
ExpAlig.exist = 0;
if ispc
    if exist([com '.exe'],'file')
        ExpAlig.exist = 1;
    end
else
    if exist(com, 'file')
        ExpAlig.exist = 1;
    end
end

confs.vars = {'flipTemplate', 'subdivLevel', 'spharmDegree', ...
        'regTemplate', 'FinalFlip', 'OutDirectory', 'others'};
confs.args = [100 1 1 100 1 200 20];
confs.flipTemplate = [SpharmMatDir DataFolder '/Ex0503/hip05_temp/a01_l_hippo_fix_surfSPHARM_ellalign.coef'];
confs.subdivLevel = 20;
confs.spharmDegree = 12;
confs.regTemplate = [SpharmMatDir DataFolder '/Ex0503/hip05_temp/a01_l_hippo_fix_surfSPHARM_ellalign.meta'];
confs.FinalFlip = 0;
confs.OutDirectory = [SpharmMatDir DataFolder '/Ex0503/hip06_coef'];
confs.others = '';
confs.command = ExpAlig.command;
confs.path = PDM_path;

if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end

if (ExpAlig.exist == 0)
    errordlg(sprintf('%s function does not exist in the PDM directory, %s',PDM_Fix.command,PDM_path));
    rmpath(Codedir);
    return;
else
    SpharmMatExpandAlignByPDM(confs, inNames);
end


%% (4) Format Conversion meta_coef2des
inFiles = dir([confs.OutDirectory '/*_ellalign.meta']); inNames={};
for i=1:length(inFiles)
    inNames{end+1} = [confs.OutDirectory '/' inFiles(i).name];
end

inFiles = dir([confs.OutDirectory '/*_ellalign.coef']); 
for i=1:length(inFiles)
    inNames{end+1} = [confs.OutDirectory '/' inFiles(i).name];
end

confs.OutDirectory = [SpharmMatDir DataFolder '/Ex0503/hip07_reg'];

if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end

SpharmMatUtilFormatConvert(confs, inNames, 'ellalign_meta_coef2reg');


%% (5) Display output objects (optional)
inFiles = dir([confs.OutDirectory '/*.mat']); inNames={};
for i=1:length(inFiles)
    inNames{end+1} = [confs.OutDirectory '/' inFiles(i).name];
end

    %Available values for Space- 'object';'param';'both'
dispConfs.Space = 'object';
    % Available values for Mesh- 'orig';'quad32';'quad64';'quad128';'quad256';'quad512';'icosa1';'icosa2'; ...
    %    'icosa3';'icosa4';'icosa5';'icosa6'
dispConfs.Mesh = 'quad64';  
    % Available values for Shape- 'solid';'mesh';'both'
dispConfs.Shade = 'both';
    % Available values for Overlay- 'none';'adc_paramap'
dispConfs.Overlay = 'adc_paramap';
    % Available values for Export- 'screen';'png';'both'
dispConfs.Export = 'png';
dispConfs.Degree = 12;
dispConfs.Template = '';

SpharmMatUtilDisplayObjs(dispConfs, inNames, CodeDir);


dispConfs.Degree = 1;
SpharmMatUtilDisplayObjs(dispConfs, inNames, CodeDir);

rmpath(CodeDir);

clear;