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
% surf_spharm.m
%
% purpose: generate surface data structure based on spharm
%
% Li Shen
% 12/29/2006 - create

function [vs, fs] = surf_spharm(fvec,dg,meshsize);

vs = []; fs = []; gc = [];

if isempty(fvec)
    return;
end

% adjust degree based on fvec
max_d = sqrt(size(fvec,1))-1;
if (dg(1)<0 | dg(2)>max_d)
    odg = dg;
    dg(1) = max(dg(1),0);
    dg(2) = min(dg(2),max_d);
    disp(sprintf('degree [%d %d] adjusted to the possible range [%d %d]', odg, dg));
else
    disp(sprintf('degree [%d %d]', dg));
end

res = meshsize;
[vs, fs] = sample_Param_Space(res);
Z = calculate_SPHARM_basis(vs, dg(2));

lb = dg(1)^2+1;
ub = (dg(2)+1)^2;

vs = real(Z(:,lb:ub)*fvec(lb:ub,:));
    
return;
