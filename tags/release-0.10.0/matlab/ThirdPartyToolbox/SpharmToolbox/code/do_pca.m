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
% do_pca.m
%
% Principle Component Analysis
%   Points are transformed by translation (zero-mean) and dimension-reduced-rotation
%
% Li Shen 
% 07/12/2002 - create

function [pca_ps, pca_b, var_amt, latent] = do_pca(points,pca_dim)

% Reduce the number of eigenvectors to be N-c (subjects minus classes).
% Then project onto the new basis, i.e. calculate the dot product of each existing
% subject vector with each of the eigenvectors.

pca_dim = min(pca_dim,size(points,2));
[pca_b, score, latent, tsquare] = princomp(points);
pca_b = pca_b(:,1:pca_dim);
pca_ps = points*pca_b;
var_amt = sum(latent(1:pca_dim))/sum(latent);
neg_ev_ind = find(latent<0);
if (~isempty(neg_ev_ind))
    disp(['Negative eigenvalues:', sprintf(' %f',latent(neg_ev_ind))]);
end 
latent = latent(1:pca_dim);

return;
