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

% Goal: 
% Input: 
%       inObjectName - File name of an input object
%       resampleFactor - Resampling factor to down- or up-sample a binary
%       object
%       saveBIM - Control flag to determine whether binary volumetric data should be
%       saved as a file for the later use. Its name will be the same as the
%       inObjectName with a different postfix('_bim.mat')
%
% Output:
%       bim - a 3D matrix of binary voxels.
%       mins - lowest coordinate values in x, y, and z axis (1 x 3)
%
% Created by Sungeun Kim (10-09-08)

function [bim, mins] = conv2bim(inObjectName, resampleFactor, saveBIM)

% Create place holders for output 
bim = [];
mins = [];

if ~ischar(inObjectName)
    disp('inObjectName should a string of input filename or list of input filenames');
    return;
elseif ~exist(inObjectName, 'file')
    disp('Input object file does not exist');
    return;
end

% Parse input object filename
[path,name,ext] = fileparts(inObjectName);

if strcmp(lower(ext), '.txt')
    token = 'TRACE';
elseif strcmp(lower(ext), '.hdr') | strcmp(lower(ext), '.nii')
    token = 'NIfTI';
else
    disp('This file format can not be processed');
    return;
end

switch (token)
    case 'TRACE'
%        [bim, mins] = zroi2bim(inObjectName, resampleFactor);  %zroi2bim
%        is incomplete and need to be further modified (10-10-08)
    case 'NIfTI'
        [bim, mins] = NIfTI2bim(inObjectName, resampleFactor);    
end

% Save binary volumetric object to a new file
if (saveBIM)
    new_name = [name '_bim'];
    save(fullfile(path, new_name), 'bim', 'mins');
end

return;