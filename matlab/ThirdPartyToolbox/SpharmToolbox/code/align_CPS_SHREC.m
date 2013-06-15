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

function [vertices, sph_verts, faces, fvec, new_name]=align_CPS_SHREC(filename, confs)
global fact;
global sgm;

gran = confs.GroupAlpha;  % gran = # of alphas to process together
res = [confs.BaseRes confs.HierarchyStep confs.HierarchyDepth confs.Top_K confs.GammaRes]; 
% base res R + step of hierarchy Hs + depth of hierarchy Hd + top N + 3rd Angle res (gammares)

switch upper(deblank(char(confs.NormalizeSize)))
    case 'YES'
        scale = 1;
    case 'NO'
        scale = 0;
end

% factorial(170) = Inf
for i=0:170 
    fact(i+1) = factorial(i);
end

utl_sgm(15);

% Load a template object
load(confs.Template);
[fvec, max_d] = fixed_fvec(fvec,confs.MaxSPHARMDegree,scale);
%tvertices = vertices; tsph_verts=sph_verts; tfvec = fvec;
tfvec = fvec;

load(filename);
[fvec, max_d] = fixed_fvec(fvec,confs.MaxSPHARMDegree,scale);
[path,name,ext] = fileparts(filename);

if ~exist('faces', 'var') | ~exist('vertices', 'var') | ~exist('sph_verts', 'var') | ~exist('fvec', 'var')
    disp('One or more of faces, vertices, spherical vertices, or SPHARM descriptor are missing');
    return;
end

rmsd_org = SPHARM_rmsd(fvec, tfvec);

if rmsd_org > 0

    % create samples in rotation space
    [alpha, beta, gamma] = utl_eas(res); % euler_angle hierarchy

    % initial alignment using ICP
    [fvec_icp, Robj_icp] = match_icp(fvec, tfvec, max_d);
    [fvec_icp, Aprm] = match_param_hie(fvec_icp,tfvec,alpha,beta,gamma,gran,res,max_d);
    [fvec_icp, Robj_icp] = match_cps(fvec_icp, tfvec, max_d);

    % assume that fvec and atlas are roughly aligned in the object space
    [fvec_cps, Robj_cps] = match_cps(fvec, tfvec, max_d);

    rmsd_icp = SPHARM_rmsd(fvec_icp, tfvec);
    rmsd_cps = SPHARM_rmsd(fvec_cps, tfvec);

    disp(sprintf('RMSD: org %0.3f, icp %0.3f, cps %0.3f',rmsd_org,rmsd_icp,rmsd_cps))

    if rmsd_org > min(rmsd_icp,rmsd_cps)    
        if rmsd_icp < rmsd_cps
            fvec = fvec_icp;
            Robj = Robj_icp;
        else
            fvec = fvec_cps;
            Robj = Robj_cps;
        end
    end
    
    % assume that fvec and atlas are roughly aligned in both object and
    % parameter spaces
    for k = 1:5
        rmsd(1) = SPHARM_rmsd(fvec, tfvec);

        % use cps to align objects together first
        [fvec_a, Robj] = match_cps(fvec, tfvec, max_d);
        rmsd(2) = SPHARM_rmsd(fvec_a, tfvec);

        if rmsd(2)<rmsd(1)
            fvec = fvec_a; % update
        end    

        % assume aligned in object space, register the parameterization
        [fvec, Aprm] = match_param_hie(fvec,tfvec,alpha,beta,gamma,gran,res,max_d);
        rmsd(3) = SPHARM_rmsd(fvec, tfvec);

        disp(sprintf('Base Res %d, Step %d, Depth %d, Top %d, Gamma %d: RMSD %0.3f => %0.3f => %0.3f',res,rmsd));

        if sum(abs(Aprm(:)))==0
            break;
        end
    end
else
    disp(sprintf('RMSD = %d: Individual (%s) is the same as template (%s)',rmsd_org,filename,confs.Template));
end

new_name = sprintf('%s/%sSHREC_reg.mat',confs.OutDirectory, name(1:end-3));
if exist(new_name,'file')
    prompt = {'Enter new filename:'};
    dlg_title = 'New File Name';
    num_lines = 1;
    def = {new_name};
    answer = inputdlg(prompt,dlg_title,num_lines,def);    
    new_name = answer{1};
end
save(new_name, 'fvec');

clear('fact','sgm');

return;
