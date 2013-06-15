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

function performTstat(confs, objs1, objs2, cpath)

% Load an atlas. If atlas name is empty, create an atlas using objs1
if exist(confs.Atlas,'file')
    load(confs.Atlas);
    atlas_fvec = fvec;
else
    conf.OutputName = 'atlas.mat';
    conf.OutDirectory = confs.OutDirectory;
    conf.SampleMesh = confs.SampleMesh;
    
    atlas_fvec = averageObjs(conf, objs1);
end

groupID=confs.GroupIDs;

% Load a sampling mesh
switch deblank(char(confs.SampleMesh))
%     case 'quad32'
%         load('ParamMeshes/R32_quad.mat');
%     case 'quad64'
%         load('ParamMeshes/R64_quad.mat');
%     case 'quad128'
%         load('ParamMeshes/R128_quad.mat');
%     case 'quad256'
%         load('ParamMeshes/R256_quad.mat');
%     case 'quad512'
%         load('ParamMeshes/R512_quad.mat');
    case 'icosa1'
        load(fullfile(cpath,'ParamMeshes/L1_icosa.mat'));
    case 'icosa2'
        load(fullfile(cpath,'ParamMeshes/L2_icosa.mat'));
    case 'icosa3'
        load(fullfile(cpath,'ParamMeshes/L3_icosa.mat'));
    case 'icosa4'
        load(fullfile(cpath,'ParamMeshes/L4_icosa.mat'));
    case 'icosa5'
        load(fullfile(cpath,'ParamMeshes/L5_icosa.mat'));
    case 'icosa6'
        load(fullfile(cpath,'ParamMeshes/L6_icosa.mat'));
end

% Reconstruct the atlas using the loaded sampling mesh
max_d = sqrt(size(atlas_fvec,1))-1;
dg = [0 max_d];
Z = calculate_SPHARM_basis(sph_verts, max_d);
lb = 1;
ub = (max_d+1)^2;
atlas_vertices = real(Z(:,lb:ub)*atlas_fvec(lb:ub,:));

% create neighborhood
nbr = get_nbr(faces,sph_verts);

% get vertex normals
vtnorm = get_vertexnormals(atlas_vertices,faces);

% Extract signal from all individuals
[defms, sgdfms, sgnrms, sgpcas, sgflds, grInfo]=extractSignals(atlas_vertices, faces, Z, max_d, vtnorm, objs1, objs2)

% Smoothe signals
FWHM = confs.Smoothing_FWHM;

for i=1:length(sgdfms)
    sgdfms{i} = performHKsmooth(sgdfms{i},atlas_vertices,faces,nbr,FWHM);
    sgnrms{i} = performHKsmooth(sgnrms{i},atlas_vertices,faces,nbr,FWHM);
    sgflds{i} = performHKsmooth(sgflds{i},atlas_vertices,faces,nbr,FWHM);
    sgpcas{i} = performHKsmooth(sgpcas{i},atlas_vertices,faces,nbr,FWHM);
end

% Perform t-test
signal = deblank(char(confs.Signal));
evar = deblank(char(confs.EqualVariance));
if strcmp(upper(evar),'YES')
    samevar = 1;
else
    samevar = 0;
end

c1sbs = find(grInfo == 1);
c2sbs = find(grInfo == 2);

switch lower(signal)
    case 'vl_defm_org'
        [tstats,pvalue] = grp_ttest(c1sbs,c2sbs,sgdfms,samevar)        
    case 'vl_defm_nrm'
        [tstats,pvalue] = grp_ttest(c1sbs,c2sbs,sgnrms,samevar)                
    case 'vl_defm_pca'
        [tstats,pvalue] = grp_ttest(c1sbs,c2sbs,sgpcas,samevar)                
    case 'vl_defm_fld'
        [tstats,pvalue] = grp_ttest(c1sbs,c2sbs,sgflds,samevar)                
end

% Save results
IDS = strrep(confs.GroupIDs,',','v');
IDS = strrep(IDS,' ','v');

new_name = [confs.OutDirectory '/' confs.OutputName];
if exist(new_name,'file')
    prompt = {'Enter new filename:'};
    dlg_title = 'New File Name';
    num_lines = 1;
    def = {new_name};
    answer = inputdlg(prompt,dlg_title,num_lines,def);    
    new_name = answer{1};
end

save(new_name, 'sph_verts', 'faces', 'atlas_vertices', 'vtnorm','grInfo','groupID','signal','tstats','pvalue','FWHM');

return;


%
% find t and p
%

function  [tstat,pvalue] = grp_ttest(c1sbs,c2sbs,sigs,samevar); % CN VS MCI

for i=1:length(c1sbs)
    id = c1sbs(i);
    ndefms1(:,:,i) = sigs{id};
end
for i=1:length(c2sbs)
    id = c2sbs(i);
    ndefms2(:,:,i) = sigs{id};
end

%  the standard deviations are unknown and not assumed equal.
d = size(ndefms1);
for i=1:d(1)
    for j=1:d(2)
        if samevar
            [h,significance,ci,stats] = ttest2(ndefms1(i,j,:),ndefms2(i,j,:),0.05);
            pvalue(i,j) = significance; tstat(i,j) = stats.tstat;
        else
            [tstat(i,j),pvalue(i,j)] = ttest2_var2(ndefms1(i,j,:),ndefms2(i,j,:));
        end
    end
end

return;

%
% two-sample t-test for equal means, variances not assumed equal.
%

function [t,p] = ttest2_var2(xs,ys)

mu1 = mean(xs); mu2 = mean(ys);
var1 = var(xs); var2 = var(ys);
n1 = length(xs); n2 = length(ys);

% t statistic
t = (mu1 - mu2)/sqrt(var1/n1 + var2/n2);
% degree of freedom
v = (var1/n1 + var2/n2)^2/((var1/n1)^2/(n1-1) + (var2/n2)^2/(n2-1));

% tail == 1
p = 1 - tcdf(t,v);
% adjust for tail == 0
p = 2 * min(p,1-p);

return;

%
% get nbr
%

function nbr = get_nbr(tri,coord)

n_points = size(coord,1);
n_tri = size(tri,1);

% compute the maximum degree of node
degree=zeros(n_points,1);
for j=1:n_tri
degree(tri(j,:))=degree(tri(j,:))+1;
end
max_degree=max(degree);

% find out the 1st neighbor nodes
nbr=zeros(n_points,max_degree);
for i_tri=1:n_tri
    for j=1:3
        cur_point = tri(i_tri,j);
        for k=1:3
            if (j ~= k)
                nbr_point= tri(i_tri,k);
                if find(nbr(cur_point,:)==nbr_point)
                    ;
                else
                    n_nbr = min(find(nbr(cur_point,:) == 0));
                    nbr(cur_point,n_nbr) = nbr_point;
                end;
            end;
        end;
    end;
end;

return;