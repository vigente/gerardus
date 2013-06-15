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

% Goal: Extract surface data (vertices and faces) from a binary object
% file (inObjectName).
% Input:
%       inObject - File name of a binary object with '_bim.mat' or
%       '_fix.mat' postfix or struct type object, consisting of 'bim' 3D matrix and
%       'mins' vector (1x3) fields.
%
%       aspRatio - Aspect ratio x:y:z to adjust coordinates depending
%       on spatial resolution in x, y, and z axis (1 x 3).
%
%       saveSurf - Control flag to determine whether surface data should be
%       saved as a file for the later use. If inObject is a string, output 
%       file name will be the same as the inObject with a different
%       postfix('_mes.mat'). Otherwise, output file name will automatically
%       be generated.
%
% Output:
%       vertices - an array of vertices (#vectices x 3).
%
%       faces - an array of faces, consisting of vertex indices (#faces x
%       3).
%
%       mins - lowest coordinate values in x, y, and z axis (1 x 3).
%
% Original code was written by Li Shen
% Modified and renamed by Sungeun Kim (10-09-08)

function [vertices, faces, mins] = extract_surface(inObject, aspRatio, saveSurf)

% Create place holders for output 
vertices = [];
faces = [];

roi = inObject

DIM = size(roi);

% Make binary image
ind = find(roi>0); roi(ind) = 1;
ind = find(roi<1); roi(ind) = 0;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Manipulate binary objects
% Note that roi will be permuted in order to make the correct north/south poles
[vertices, faces] =  gen_surf_data(roi)
vertices = vertices + mins(ones(1,size(vertices,1)),:);

% Adjust aspect ratio
vertices = vertices.*asratio(ones(1,size(vertices,1)),:);

% Save surface object to a new file
if (saveSurf)
    if ischar(inObject)
        new_name = [name(1:end-3) 'obj'];
        save(fullfile(path, new_name), 'vertices', 'faces', 'mins');
    else
        name = genvarname('default_surface', who);
        new_name = [name '_mes'];
        save(fullfile(pwd, new_name), 'vertices', 'faces', 'mins');
    end
end

return;
