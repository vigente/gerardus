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

%
% create spherical harmonic descriptor
%

function [fvec, deg, Z, new_name] = create_SPHARM_des_LSF(vertices, faces, sph_verts ,maxDeg, filename, OutDirectory)

if isempty(vertices) | isempty(sph_verts)
    if ~isempty(filename)
        load(filename);
    else
        disp('There is no useful information');
        return;
    end
end

[pa, na, ex] = fileparts(filename);
new_name = '';

if ~exist('vertices', 'var') | ~exist('sph_verts', 'var')
    disp('One or more of vertices, and spherical vertices are missing');
    return;
end

vertnum = size(sph_verts,1);

max_d = maxDeg;
% Note that degree 'd' we want to use depends on the vertnum 
% The total number of unknowns is (d+1)*(d+1)
% The total number of equations is vertnum
% We want equ_num >= unk_num
deg = max(1, floor(sqrt(vertnum)*1/2));
deg = min(deg, max_d);
disp(sprintf('Use spharm up to %d degree (vec_len=%d).',deg,(deg+1)^2));

Z = calculate_SPHARM_basis(sph_verts, deg); 

[x,y] = size(Z);
disp(sprintf('Least square for %d equations and %d unknowns',x,y));

% Least square fitting
% fvec = Z\vertices;   %This does not work as it is expected to work in certain environment
for i=1:size(vertices,2)
    fvec(:,i) = Z\vertices(:,i);
end

if ~isempty(filename)
    if ~isempty(OutDirectory)
        new_name = sprintf('%s/%s_LSF_des.mat', OutDirectory, na(1:end-4));
    else
        new_name = sprintf('%s/%s_LSF_des.mat', pa, na(1:end-4));
    end
    if exist(new_name, 'file')
        prompt = {'Enter new filename:'};
        dlg_title = 'New File Name';
        num_lines = 1;
        def = {new_name};
        answer = inputdlg(prompt,dlg_title,num_lines,def);    
        new_name = answer{1};
    end
    save(new_name, 'vertices', 'faces', 'sph_verts', 'fvec');
end

return;
