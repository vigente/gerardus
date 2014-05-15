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
%
%

function write_metaNcoef(filename, confs)

load(filename);

if ~exist('faces', 'var') | ~exist('vertices', 'var') | ~exist('fvec', 'var')
    disp('One or more of faces, vertices, and SPHARM coefficients are missing');
    return;
end

[fvec] = read_coef(coef);

[pa, na, ex]=fileparts(filename);

type = size(faces,2);
filename1 = sprintf('%s/%s_surfSPHARM.meta', confs.OutDirectory, na(1:end-4));
filename2 = sprintf('%s/%s_surfSPHARM.coef', confs.OutDirectory, na(1:end-4));

if exist(filename1,'file') & exist(filename2,'file')
    prompt = {'Enter new filenames:'};
    dlg_title = 'New File Names';
    num_lines = 2;
    def = {filename1; filename2};
    answer = inputdlg(prompt,dlg_title,num_lines,def);    
    filename1 = answer{1};
    filename2 = answer{2};
end

% Write a meta file
if type == 3
    disp('This function does not support any objects in the general triangular mesh now');
elseif type == 4
    write_meta_quad(filename1, vertices, faces);
%    write_meta_quad(filename2, sph_verts, faces);
end

% Write a coef file
write_coef(filename2, fvec);

return;