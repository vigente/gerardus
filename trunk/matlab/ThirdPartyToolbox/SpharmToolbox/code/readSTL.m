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

function [vertices, faces] = readSTL(filename)

[fs, vs, cout] = stlreadASCII(filename);

[v,I,J] = unique(vs,'rows');

vertices = vs(sort(I),:);
faces = zeros(size(fs));
for i = 1:size(fs,2)
    vf = vs(fs(:,i), :);
    [tf, lc] = ismember(vf,vertices, 'rows');
    faces(:,i) = lc;
end

faces = unique(faces,'rows');
dif1 = faces(:,1)-faces(:,2); dif2=faces(:,2)-faces(:,3); dif3=faces(:,3)-faces(:,1);
indDif1 = find(dif1 == 0);
indDif2 = find(dif2 == 0);
indDif3 = find(dif3 == 0);
indDif = union(indDif1, indDif2, indDif3);

ufIDX = setdiff([1:size(faces,1)], indDif);
faces = faces(ufIDX,:);

return;