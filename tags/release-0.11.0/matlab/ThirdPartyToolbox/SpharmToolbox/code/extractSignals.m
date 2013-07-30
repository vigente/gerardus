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

function [defms, vldfms, vlnrms, vlpcas, vlflds, grInfo] = extractSignals(atlas, faces, Z, max_d, vnorm, objs1, objs2)

lb = 1;
ub = (max_d+1)^2;

numSbj1 = length(objs1);
numSbj2 = length(objs2);

defms=cell(numSbj1+numSbj2,1);  
vldfms=cell(numSbj1+numSbj2,1);  vlnrms=cell(numSbj1+numSbj2,1);  
vlpcas=cell(numSbj1+numSbj2,1);  vlflds=cell(numSbj1+numSbj2,1);
grInfo = ones(numSbj1+numSbj2,1);

% Group 1

for i=1:numSbj1
    file = objs1{i};
    [pa,na,ex]=fileparts(file);

    load(file);
    % Reconstruction an individual
    vertices = real(Z(:,lb:ub)*fvec(lb:ub,:));

    % Calculate deformation fields
    defm = vertices - atlas;
    
    vldefm = sqrt(sum(defm.^2,2));
    vlnorm = sum(defm.*vnorm,2);

    negSGN = find(vlnorm<0);
    vldefm(negSGN) = vldefm(negSGN) * -1;

    defms{i} = defm;
    vldfms{i} = vldefm;
    vlnrms{i} = vlnorm;
    
    grInfo(i) = 1;
end

% Group 2

for i=1:numSbj2
    file = objs2{i};
    [pa,na,ex]=fileparts(file);

    load(file);
    % Reconstruction an individual
    vertices = real(Z(:,lb:ub)*fvec(lb:ub,:));

    % Calculate deformation fields
    defm = vertices - atlas;
    
    vldefm = sqrt(sum(defm.^2,2));
    vlnorm = sum(defm.*vnorm,2);

    negSGN = find(vlnorm<0);
    vldefm(negSGN) = vldefm(negSGN) * -1;

    defms{i+numSbj1} = defm;
    vldfms{i+numSbj1} = vldefm;
    vlnrms{i+numSbj1} = vlnorm;
    
    grInfo(i+numSbj1) = 2;
end

% calculate PCA and FLD direction
c1sbs = find(grInfo == 1);
c2sbs = find(grInfo == 2);
[vpca, vfld] = get_vprin(c1sbs,c2sbs,atlas,faces,defms,vnorm);

% calculate sgpca and sgfld
for i = 1:numSbj1;
	vlpca = sum(defms{i}.*vpca,2);
	vlfld = sum(defms{i}.*vfld,2);
    vlpcas{i}=vlpca;  
    vlflds{i}=vlfld;
end

for i = 1:numSbj2;
	vlpca = sum(defms{i+numSbj1}.*vpca,2);
	vlfld = sum(defms{i+numSbj1}.*vfld,2);
    vlpcas{i+numSbj1}=vlpca;  
    vlflds{i+numSbj1}=vlfld;
end

return;


%
% template: mean of CNs
% defms: deformation fields from the template (CN mean) 
% norms: difference along surface normal directions
%

function [vpca, vfld] = get_vprin(c1sbs,c2sbs,template,fs,defms,vnorm)

% consider only outward PCA direction
[vpca, vfld] = get_directions(c1sbs,c2sbs,defms);

% fix vpca and vfld to get outward direction
idx = find(sum(vpca.*vnorm,2)<0);
vpca(idx,:) = vpca(idx,:)*-1;
vnlen = sqrt(sum(vpca.^2,2));
vpca = vpca./vnlen(:,[1 1 1]);

idx = find(sum(vfld.*vnorm,2)<0);
vfld(idx,:) = vfld(idx,:)*-1;
vnlen = sqrt(sum(vfld.^2,2));
vfld = vfld./vnlen(:,[1 1 1]);

return;

%
% get PCA and FLD directions (length one)
%

function [vpca, vfld] = get_directions(c1sbs,c2sbs,defms)

n1 = length(c1sbs);
for i = 1:n1
    c1defms(:,:,i) = defms{c1sbs(i)};
end
n2 = length(c2sbs);
for i = 1:n2
    c2defms(:,:,i) = defms{c2sbs(i)};
end

d = size(c1defms);
for i = 1:d(1)
    vs1 = squeeze(c1defms(i,:,:))';
    vs2 = squeeze(c2defms(i,:,:))';

    vs = [vs1;vs2];
    [pca_ps, pca_b, var_amt, latent] = do_pca(vs,3);
    
    Labels = [ones(n1,1); ones(n2,1)*2];
    [FLD_basis, FLD_vals] = FLD(vs,Labels);
    
    PCA_basis = pca_b(:,1); PCA_basis = PCA_basis/norm(PCA_basis);
    FLD_basis = FLD_basis/norm(FLD_basis);

    vpca(i,:) = PCA_basis';
    vfld(i,:) = FLD_basis';
    
end

return;