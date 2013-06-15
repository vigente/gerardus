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

function write_procalign(filename, confs)

load(filename);

if ~exist('fvec', 'var') 
    disp('SPHARM descriptor is missing');
    return;
end

[pa, na, ex, ve]=fileparts(filename);

strMesh = confs.SampleMesh;
switch upper(strMesh)
    case 'quad32'
        res = 32;
    case 'quad64'
        res = 64;
    case 'quad128'
        res = 128;    
    case 'quad256'
        res = 256;        
    case 'quad512'
        res = 512;        
    case 'icosa1'
        res = -1;        
    case 'icosa2'
        res = -2;        
    case 'icosa3'
        res = -3;        
    case 'icosa4'
        res = -4;        
    case 'icosa5'
        res = -5;        
    case 'icosa6'
        res = -6;        
end

% Resample (reconstruct individuals from SPHARM descriptor
[sph_verts, faces] = sample_Param_Space(res);
maxDeg = sqrt(size(fvec,1))-1;
Z = calculate_SPHARM_basis(sph_verts, maxDeg);
lb = 1;   ub = maxDeg;
vertices = real(Z(:,lb:ub)*fvec(lb:ub,:));

type = size(faces,2);
filename1 = sprintf('%s/%s_procalign.meta', confs.OutDirectory, na(1:end-4));

if exist(filename1,'file')
    prompt = {'Enter new filename:'};
    dlg_title = 'New File Name';
    num_lines = 1;
    def = {filename1};
    answer = inputdlg(prompt,dlg_title,num_lines,def);    
    filename1 = answer{1};
end

if type == 3
    write_meta_tri(filename1, vertices, faces);
elseif type == 4
    write_meta_quad(filename1, vertices, faces);
end

return;
