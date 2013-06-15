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
% write_gipl.m
%
% convert *.mat file to *.gipl file with scales = [1 1 1]
% 
% Wan Jing
% 11/06/2008 - create
%============================================

function write_gipl(fname, confs) 
% fname is the name of the input mat file. 
% If not specifying a filename, a dialogue window will be activated for choosing a
% mat file. 

if ~exist(fname, 'file')
    disp('Input object doesnot exist');
    return;
end

load(fname);

scales = vxsize;
if max(size(origin)==3)
    origin(4)= 0;
end

[path, name, ext,ver] = fileparts(fname);

postfix = name(end-2:end);
if strcmp(postfix, 'bim') | strcmp(postfix, 'fix')
    new_name = [confs.OutDirectory '/' name(1:end-4) '.gipl'];
else
    new_name = [confs.OutDirectory '/' name '.gipl'];
end
if exist(new_name,'file')
    prompt = {'Enter new filename:'};
    dlg_title = 'New File Name';
    num_lines = 1;
    def = {new_name};
    answer = inputdlg(prompt,dlg_title,num_lines,def);    
    new_name = answer{1};
end
            
oldsize = size(bim);
bim_new = zeros(oldsize(1)+2, oldsize(2)+2, oldsize(3)+2);
bim_new(2:end-1, 2:end-1, 2:end-1) = bim;
origin = origin-1;
% generate a gipl file
gipl_write_volume(bim_new,new_name,scales,origin);

return;
