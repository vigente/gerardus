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
% FLD.m
%
% Goal: 
%   1. Calculate the optimal FLD projection
%   2. For two class problem
%
% Li Shen 
% 11/13/2002 - create

function [FLD_basis, FLD_vals] = FLD(Samples,Labels)

d = size(Samples);
if (d(2)>d(1)-2)
    disp('---------------------------------');
    disp(sprintf('N=%d points, M=%d dims (IGNORE the last %d dims): a nonsigular Sb requires M<=N-c (c=2)',d,d(2)-d(1)+2));
    disp('---------------------------------');
    Samples = Samples(:,1:d(1)-2);
end

FLD_basis = [];

Sb = get_Sb(Samples,Labels);
Sw = get_Sw(Samples,Labels);

[V,D] = eig(Sb,Sw);

[eigval,ind] = sort(diag(D));

FLD_basis = V(:,ind(end));

FLD_vals = Samples*FLD_basis;

return;

%
% calculate between class scatter matrix Sb
%

function Sb = get_Sb(Samples,Labels)

c1_ind = find(Labels==1); c2_ind = find(Labels==2);
m = mean(Samples,1); m1_m = mean(Samples(c1_ind,:),1) - m ; m2_m = mean(Samples(c2_ind,:),1) - m;
Sb = length(c1_ind)*(m1_m'*m1_m) + length(c2_ind)*(m2_m'*m2_m);

return;

%
% calculate within class scatter matrix Sw
%

function Sw = get_Sw(Samples,Labels)

c1_ind = find(Labels==1); c2_ind = find(Labels==2);

% need to have at least 2 points in each class
Sw = cov(Samples(c1_ind,:))*(length(c1_ind)-1) + cov(Samples(c2_ind,:))*(length(c2_ind)-1);

return;

%
% calculate within class scatter matrix Sw (directly according to the definition)
%

function Sw = get_Sw_v2(Samples,Labels)

c1_ind = find(Labels==1); c2_ind = find(Labels==2);

% to verify the correctness of Sw
Sw = zeros(size(Samples,2));
m1 = mean(Samples(c1_ind,:),1);
xs1 = Samples(c1_ind,:) - m1(ones(1,length(c1_ind)),:);
for i = 1:size(xs1,1)
    Sw = Sw + xs1(i,:)'*xs1(i,:);
end

m2 = mean(Samples(c2_ind,:),1);
xs2 = Samples(c2_ind,:) - m2(ones(1,length(c2_ind)),:);
for i = 1:size(xs2,1)
    Sw = Sw + xs2(i,:)'*xs2(i,:);
end

return;