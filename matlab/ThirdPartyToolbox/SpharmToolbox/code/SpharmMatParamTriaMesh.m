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

function [vertices, faces, sph_verts, new_name] = parameterizeCALD(file, confs)

load(file);
[path, name, ext] = fileparts(file);

% First, if obj == bim, then extract surface
if exist('bim','var')
    % Make binary image
    ind = find(bim>0); bim(ind) = 1;
    ind = find(bim<1); bim(ind) = 0;    
    
    [vertices, faces] =  gen_surf_data(bim,origin,vxsize);
    
    postfix = name(end-2:end);
    if strcmp(postfix,'bim') | strcmp(postfix,'fix')
        name2 = [path '/' name(1:end-3) 'obj.mat'];
    else
        name2 = [path '/' name '_obj.mat'];
    end
else
    name2 = file;
end

% Second, perform the initial parameterization
[sph_verts, name3] = initParamCALD(vertices, faces, name2, confs);

% Third, smooth the parameterization
[vertices, faces, sph_verts, new_name] = smootheCALD(vertices, faces, sph_verts, name3, confs);

return;