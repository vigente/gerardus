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
% calculate area scaling ratio between the object mesh and the parameter mesh
%

function asr = calc_asr(vertices, sph_verts, faces)

% convert quadralaterals to triangles
if size(faces,2)==4
    disp('Please use triangular meshes!');
    asr.maxcon = nan;
    asr.maxexp = nan;
    asr.objstc = nan; 
    asr.prmstc = nan;
    return;
%     faces = [faces(:,1:3); faces(:,[3 4 1])];
end

% calculate relative areas of triangles on object surface net  
[obj_area, count] = calc_triangle_areas(vertices, faces, 'object', 'triangle');
% calculate relative areas of spherical triangles in parameter space
[par_area, count] = calc_triangle_areas(sph_verts, faces, 'parameter', 'triangle');

val = par_area./obj_area; stretch = max(val,1./val);

asr.maxcon = max(1./val);
asr.maxexp = max(val);
asr.objstc = sum(stretch.*obj_area); 
asr.prmstc = sum(stretch.*par_area); % use this as area distortion cost

return;
