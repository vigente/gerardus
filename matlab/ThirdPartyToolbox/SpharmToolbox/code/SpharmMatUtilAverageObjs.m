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

function [fvec] = SpharmMatUtilAverageObjs(confs, objs)   
global fact;
global alpha;
global beta;
global gamma;
global gran;
global res;
global dg;

atlasName = confs.OutputName;
outDir = confs.OutDirectory;
if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end
if (~exist([confs.OutDirectory '/Logs'],'dir'))
    mkdir([confs.OutDirectory '/Logs']);
end

diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' atlasName '_avgObjs.log']));

numSbj = size(objs,2);

% factorial(170) = Inf
for i=0:170 
    fact(i+1) = factorial(i);
end

utl_sgm(15);

% inform = spm_input('gran scale', '+0', 's', '100 0'); 
inform = '100 0'; 
info = str2num(inform);
gran = info(1); % gran = # of alphas to process together
scale = info(2); % make rmsd scaling invariant
inform = '1 1 3 1 2';
res = str2num(inform); % base res R + step of hierarchy Hs + depth of hierarchy Hd + top N + 3rd Angle res (gammares)

% create samples in rotation space
[alpha, beta, gamma] = utl_eas(res); % euler_angle hierarchy

% collect the data
for i = 1:numSbj
    file = objs{i};
    [pa,na,ex,ve]=fileparts(file);
    newfile{i} = fullfile(confs.OutDirectory,[na(1:end-3) 'als' ex]);

    load(file);
    if ~exist('fvec','var')
        disp('Current version requires SPHARM descriptors');
        break;
    end

    max_degree = sqrt(size(fvec,1))-1;
    [fvec,max_d] = fixed_fvec(fvec,max_degree,0);
    fvecs(:,:,i) = fvec;
end


% initialize the atlas
atlas = fvecs(:,:,1);
dg = [0 max_d]; 
     
% repeat until convergence
k = 1;
while 1
    for i=1:numSbj
        fvec = fvecs(:,:,i);
        [fvec,rmsd(i,:)] = shrec_kth_itr(fvec,atlas,max_d,k);
        fvecs(:,:,i) = fvec;
        save(newfile{i}, 'fvec','k');
	end

    prev = atlas;
    atlas = mean(fvecs,3); 

    d = SPHARM_rmsd(prev, atlas);
    
    if abs(d)<0.01
        disp(sprintf('Successful: k=%d, d=%f',k,d));
        break;
    end
    k = k+1;
end

fvec = atlas;

new_name = [outDir '/' atlasName];
if exist(new_name,'file')
    prompt = {'Enter new filename:'};
    dlg_title = 'New File Name';
    num_lines = 1;
    def = {new_name};
    answer = inputdlg(prompt,dlg_title,num_lines,def);    
    new_name = answer{1};
end

save(new_name, 'fvec');
diary('off');

clear('fact','alpha','beta','gamma','gran','res','dg');
return;