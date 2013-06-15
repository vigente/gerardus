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
% write_coef.m
% function for writing *.coef file 

% INPUT: filename, which is the name of the input coef file. 
% If not specifying a filename, a dialogue window will be activated for choosing a
% gipl file. 
%
% Sungeun
%============================================

function write_coef(filename, fvec)

fid = fopen(filename,'w');
if(fid<0)
    fprintf('could not create file %s\n',filename);
    return;
end


Count = size(fvec,1);
fprintf(fid,'{ %d,',Count);

% write the coefficients to output file
idx = 0;
md = sqrt(Count)-1;
for L = 0:md
    for M=-L:L
        idx = idx+1;
        if idx < Count
            switch sign(M)
                case -1

                case 0
                    fprintf(fid, '{%lf, %lf, %lf},\n', real(fvec(idx,1)),real(fvec(idx,2)),real(fvec(idx,3)));
                case 1
                    fprintf(fid, '{%lf, %lf, %lf},\n', real(fvec(idx,1)),real(fvec(idx,2)),real(fvec(idx,3)));
                    fprintf(fid, '{%lf, %lf, %lf},\n', imag(fvec(idx,1)),imag(fvec(idx,2)),imag(fvec(idx,3)));
            end
        else
                fprintf(fid, '{%lf, %lf, %lf},\n', real(fvec(idx,1)),real(fvec(idx,2)),real(fvec(idx,3)));
                fprintf(fid, '{%lf, %lf, %lf}}', imag(fvec(idx,1)),imag(fvec(idx,2)),imag(fvec(idx,3)));
        end
    end
end

fclose(fid);

return