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

disp('Ex0701: Import objects')

SpharmMatDir = '../';
CodeDir = [SpharmMatDir 'code'];
addpath(CodeDir);

DataFolder = 'data';

%% (1) Input data directory
InDataDir = [SpharmMatDir DataFolder '/Ex0701/hip00_org'];


%% (2) List input file names
inFiles = dir([InDataDir '/*_surf.m']); inNames1={};
for i=1:length(inFiles)
    inNames1{end+1} = [InDataDir '/' inFiles(i).name];
end

inFiles = dir([InDataDir '/*_bim.hdr']); inNames2={};
for i=1:length(inFiles)
    inNames2{end+1} = [InDataDir '/' inFiles(i).name];
end

inFiles = dir([InDataDir '/*_surf.stl']); inNames3={};
for i=1:length(inFiles)
    inNames3{end+1} = [InDataDir '/' inFiles(i).name];
end

%% (3) Import Objects
confs.ResampleFactor=1; 
confs.OutDirectory = [SpharmMatDir DataFolder '/Ex0701/hip01_res'];

if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end

SpharmMatUtilImportObjs(confs, inNames1);
SpharmMatUtilImportObjs(confs, inNames2);
SpharmMatUtilImportObjs(confs, inNames3);

%% (4) Display output objects (optional)
    %Available values for Space- 'object';'param';'both'
dispConfs.Space = 'object';
    % Available values for Mesh- 'orig';'quad32';'quad64';'quad128';'quad256';'quad512';'icosa1';'icosa2'; ...
    %    'icosa3';'icosa4';'icosa5';'icosa6'
dispConfs.Mesh = 'orig';  
    % Available values for Shape- 'solid';'mesh';'both'
dispConfs.Shade = 'both';
    % Available values for Overlay- 'none';'adc_paramap'
dispConfs.Overlay = 'none';
    % Available values for Export- 'screen';'png';'both'
dispConfs.Export = 'png';
dispConfs.Degree = [];
dispConfs.Template = '';

inFiles = dir([confs.OutDirectory '/*.mat']); inNames={};
for i=1:length(inFiles)
    inNames{end+1} = [confs.OutDirectory '/' inFiles(i).name];
end

SpharmMatUtilDisplayObjs(dispConfs, inNames, CodeDir);


rmpath(CodeDir);

clear;