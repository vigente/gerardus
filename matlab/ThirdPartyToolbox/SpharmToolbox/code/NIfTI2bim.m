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

% Goal: Convert NIfTI image to bim format
% Input:
%       inObjectName - File name of an input object with .nii or .hdr
%       extension
%
%       resampleFactor - Resampling factor to down- or up-sample a binary
%       object
%
% Output:
%       bim - a 3D matrix of binary voxels.
%       mins - lowest coordinate values in x, y, and z axis (1 x 3)
%
% 03/31/2001 - created by Li Shen
% Modified by Sungeun Kim (10-09-08)

function [bim, mins] = NIfTI2bim(inObjectName, resampleFactor)

ni = load_nii(inObjectName);
roi = ni.img;
clear ni;

DIM = size(roi);

if resampleFactor ~= 1   % Upsampling or downsampling
    roi = imresize(roi, resampleFactor, 'nearest');

    space = 1/resampleFactor;  % Interpolate binary volume along z axis
    if space < 1        % Up-sampling
        [xInt, yInt, zInt] = meshgrid(1:size(roi,1), 1:size(roi,2), 1:space:size(roi,3)+space);
    elseif space > 1    % Down-sampling
        [xInt, yInt, zInt] = meshgrid(1:size(roi,1), 1:size(roi,2), 1:space:size(roi,3));
    end
    bim = interp3(roi, xInt, yInt, zInt,'nearest');
else
    bim = roi;
end
    
% Make binary image
ind = find(bim>0); bim(ind) = 1;
ind = find(bim<1); bim(ind) = 0;  

mins = [0 0 0];

return;