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
% read_meta.m
%
% read *.meta file and visualize it
% 
% Wan Jing
% 11/06/2008 - create
%============================================
function  [Points,Cells2] = read_meta(filename)
% filename is the name of the input meta file. 
% If not specifying a filename, a dialogue window will be activated for choosing a
% meta file. 

fid = fopen(filename,'r');
if(fid<0)
    fprintf('could not open file %s\n',filename);
    return
end

% read the head of meta file and save the information to Info
for i = 1:20
    Line = fgetl(fid);
    for j = 1:length(Line)
        Info(i, j) = Line(j);
    end
    str = [Line 'abcdefgh'];
    %disp(sprintf('%d: %s',i,str))
    if strcmp(str(1:8),'PointDim')==1
        nline = i;
        break;
    end
end

% npoints store the number of the points readed from meta file
npoints = fscanf(fid, 'NPoints = %d\n', 1);

Line = fgetl(fid);

% read the points and save them to the array named Points. 
% each row contains x-,y-,z- coordinate value of a point.
Points = fscanf(fid, '%*d %f %f %f \n', [3, npoints]);
Points = Points';

%store head information to info
Line = fgetl(fid);
for j = 1:length(Line)
    Info(nline+1,j) = Line(j);
end

if strcmp(Line,'CellType = TRI')==1
    %store the number of cells to ncells
    ncells = fscanf(fid, 'NCells = %d\n', 1);

    %print info, npoints, ncells to screen.
%    Info
    fprintf('The number of points is:  %d\n',npoints);
    fprintf('The number of cells is:  %d\n',ncells);

    Line = fgetl(fid);

    Cells = fscanf(fid, '%*d %d %d %d \n', [3, ncells]);
    Cells = Cells';

    Cells2 = Cells+1;

    %if want to show quadrilaterals instead of triangle, activate the
    %following line
     %Cells2 = [Cells2(1:2:end,:) Cells2(2:2:end,2)];
    
elseif strcmp(Line,'CellType = QUAD')==1 
     ncells = fscanf(fid, 'NCells = %d\n', 1);
%     Info
     fprintf('The number of points is:  %d\n',npoints);
     fprintf('The number of cells is:  %d\n',ncells);
     
     Line = fgetl(fid);
     
     Cells = fscanf(fid, '%*d %d %d %d %d \n', [4, ncells]);
     Cells = Cells';

     Cells2 = Cells+1;
else
    fprintf('The cell type is not supported \n'); 
    fclose('all');
    return
end

%plot the 3-D mesh object.
%figure, patch_mesh(Points,Cells2);
fclose('all');
return


