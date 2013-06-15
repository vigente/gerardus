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
% read_coef.m
% function for reading *.coef file and saving the coefficients to *_des.mat
% file

% INPUT: filename, which is the name of the input coef file. 
% If not specifying a filename, a dialogue window will be activated for choosing a
% gipl file. 
%
% Wan Jing
% 11/06/2008 - create
%============================================

function fvec = read_coef(filename)

fid = fopen(filename,'r');
if(fid<0)
    fprintf('could not open file %s\n',filename);
    return
end


Count = fscanf(fid, '{ %d,', 1);
if(isempty(Count))
    fprintf('The type of this file is not coef. \n');
    return;
end

% read the coefficients to C
C = fscanf(fid, '{%lf, %lf, %lf},\n', [3 inf]);
C = C';
fprintf('The number of coefficients are: %d .\n ',length(C));

if (Count ~= length(C))
    fprintf('The file has error about the number of coefficients');
end

fvec = [];
i = 1;

md = sqrt(length(C))-1;
for L = 0:md
    for M=-L:L
%         disp(sprintf('i=%d, L=%d, M=%d',i,L,M));
        switch sign(M)
            case -1
                fvec(end+1,:) = [0 0 0];
            case 0
                fvec(end+1,:) = C(i,:); i=i+1;
            case 1
                fvec(end+1,:) = complex(C(i,:),C(i+1,:)); i=i+2;
        end
    end
end

fclose('all');

return

    