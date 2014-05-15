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

function scaleObjs(confs, objs)

numSbj = size(objs,2);
outDir = confs.OutDirectory;

if ~exist(outDir,'dir')
    mkdir(outDir);
end;
if (~exist([confs.OutDirectory '/Logs'],'dir'))
    mkdir([confs.OutDirectory '/Logs']);
end

if ~exist(confs.ScalingFactor,'file')
    disp('Multiplicative scaling factor information is missing');
    return;
end

[sfactors, objNames] = loadSF(confs.ScalingFactor);

for i = 1:numSbj
    file = objs{i};
    [pa,na,ex]=fileparts(file);

    diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' na '_scaleObjs.log']));
    
    sIDX = find(strcmpi([na ex], objNames)==1);
    
    if ~isempty(sIDX)
        sf = (sfactors(sIDX(1)))^(1/3);
        load(file);
        
        if ~exist('fvec','var') & ~exist('vertices','var')
            disp(sprintf('Ignore %s',file));
            continue;
        end
                
        if exist('fvec','var')
            fvec = fvec*sf;
        end
        
        if exist('vertices','var')
            vertices = vertices*sf;
        end
        
        new_name = [outDir '/' na(1:end-3) 'scl_' na(end-2:end) ex];
        if exist(new_name,'file')
            prompt = {'Enter new filename:'};
            dlg_title = 'New File Name';
            num_lines = 1;
            def = {new_name};
            answer = inputdlg(prompt,dlg_title,num_lines,def);    
            new_name = answer{1};
        end

%         *_obj.mat: vertices, faces
%         *_ini.mat, *_smo.mat: vertices, sph_verts
%         *_des.mat, *_prm.mat, *_reg.mat: fvec, vertices, sph_verts
%         *_reg.mat: fvec without vertices
        
        if ~exist('vertices','var') % *_reg.mat (without vertices)
            save(new_name,'sf', 'fvec'); 
        elseif ~exist('sph_verts','var') % *_obj.mat
            save(new_name,'sf', 'vertices','faces'); 
        elseif ~exist('fvec','var') % *_ini.mat, *_smo.mat
            save(new_name,'sf', 'vertices','sph_verts','faces');
        else % *_des.mat, *_prm.mat, *_reg.mat
            save(new_name,'sf', 'fvec','vertices','sph_verts','faces');
        end
    else
        disp(sprintf('Scaling factor does not exist for %s%s\n',na,ex));
    end
    
    diary('off');
end    
    
return;

% It is assumed about a SF file that
% (1) it is a plain text files in Comma delimited format without any header
% (2) each row in a file has a filename with a multiplicative scaling
% factor for one object. For example, name1,0.9

function [sf, fnames] = loadSF(fname)

fin = fopen(fname,'rt');
sf = []; fnames = {};
while 1
   tline = fgetl(fin);
   if ~ischar(tline), break; end;
   pos = strfind(tline,',');
   fnames{end+1} = strtrim(tline(1:pos(1)-1));
   sf(end+1) = str2num(deblank(tline(pos(1)+1:end)));
end

fclose(fin);

return;