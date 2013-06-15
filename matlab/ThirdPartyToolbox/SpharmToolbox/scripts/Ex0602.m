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

disp('Ex0602: Perform PCA')

SpharmMatDir = '../';
CodeDir = [SpharmMatDir 'code'];
addpath(CodeDir);

DataFolder = 'data';

%% (1) Input data directory
InDataDir = [SpharmMatDir DataFolder '/Ex0602/hip06_reg'];


%% (2) List input file names
inFiles = dir([InDataDir '/*_reg.mat']); inNames={};
for i=1:length(inFiles)
    inNames{end+1} = [InDataDir '/' inFiles(i).name];
end


%% (3) Display input objects (optional)
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
dispConfs.Degree = [];
dispConfs.Template = '';

SpharmMatUtilDisplayObjs(dispConfs, inNames, CodeDir);


%% (4) Perform PCA
confs.GroupID = 'All';
confs.OutputName = 'PCA_stat.mat';
confs.OutDirectory = [SpharmMatDir DataFolder '/Ex0602/hip07_stat'];

if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end

SpharmMatStatAnalysis(confs, inNames, {}, 'PCA', CodeDir);        


%% (5) Display output objects (optional)
dispConfs.Threshold_p_value = 0.05;
    % Available values for Overlay- 'p-value';'t-map'
dispConfs.Overlay = 'p-value';
    % Available values for Colormap-
    % 'jet';'hot';'summer';'cool';'autumn';'winter';'spring'
dispConfs.Colormap = 'jet';

dispConfs.Level = 4;
dispConfs.Sigma = 3;
    % Available values for Colormap- 'quad32';'quad64';'quad128';'quad256';'quad512';'icosa1';'icosa2'; ...
    % 'icosa3';'icosa4';'icosa5';'icosa6'
dispConfs.Mesh = 'quad32';
dispConfs.MaxSPHARMDegree = 15;

inFile = dir([confs.OutDirectory '/' confs.OutputName]); inName={};
for i=1:length(inFile)
    inName{end+1} = [confs.OutDirectory '/' inFile(i).name];
end

h = figure('Name', 'Render Statistical Results','NumberTitle', 'off');
cameratoolbar(h, 'Show');
SpharmMatDisplayStat(dispConfs, inName, 'res_PCA', h, CodeDir);        

rmpath(CodeDir);

clear;