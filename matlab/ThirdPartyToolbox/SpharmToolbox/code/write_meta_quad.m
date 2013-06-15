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

% ============================================
% write_meta_Quad.m
%
% write *.meta file in Quadrilateral
% 
% Wan Jing
% 11/08/2008 - create
%============================================
function write_meta_Quad(filename, vertices, faces)

fid = fopen(filename, 'w');
fprintf(fid, 'ObjectType = Scene\n');
fprintf(fid, 'NDims = 3\n');
fprintf(fid, 'NObjects = 1\n');
fprintf(fid, 'ObjectType = Mesh\n');
fprintf(fid, 'NDims = 3\n');
fprintf(fid, 'BinaryData = False\n');
fprintf(fid, 'TransformMatrix = 1 0 0 0 1 0 0 0 1\n');
fprintf(fid, 'Offset = 0 0 0\n');
fprintf(fid, 'CenterOfRotation = 0 0 0\n');
fprintf(fid, 'ElementSpacing = 1 1 1\n');
fprintf(fid, 'PointType = MET_FLOAT\n');
fprintf(fid, 'PointDataType = MET_FLOAT\n');
fprintf(fid, 'CellDataType = MET_FLOAT\n');
fprintf(fid, 'NCellTypes = 1\n');
fprintf(fid, 'PointDim = ID x y ...\n');
fprintf(fid, 'NPoints = %d\n', length(vertices));
fprintf(fid, 'Points = \n');
for i = 0:(length(vertices)-1)
    fprintf(fid,'%d %g %g %g \n', i, vertices(i+1,1), vertices(i+1,2), vertices(i+1, 3));
end
fprintf(fid, 'CellType = QUAD\n');
fprintf(fid, 'NCells = %d\n',length(faces));
fprintf(fid, 'Cells = \n'); 
for i = 0:(length(faces)-1)
    fprintf(fid, '%d %d %d %d %d \n', i, (faces(i+1,1)-1), (faces(i+1,2)-1), (faces(i+1,3)-1), (faces(i+1,4)-1));
end
fclose(fid);

return;