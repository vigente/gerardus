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

%
% use cps to align object together first
%

function [rvec, M] = match_cps(fvec,atlas,max_d)

dg = [0 max_d]; meshsize = -3;
[X, fs] = surf_spharm(atlas,dg,meshsize);
[P, fs] = surf_spharm(fvec,dg,meshsize);

% align P to X
[P,M] = align_cps(P,X);

rvec = (M(1:3,1:3)*fvec')';
rvec(1,:) = fvec(1,:) + M(1:3,4)'*2*sqrt(pi); % Y00 = 1/(2*sqrt(pi))

rmsd(1) = SPHARM_rmsd(fvec, atlas);
rmsd(2) = SPHARM_rmsd(rvec, atlas);

if rmsd(1)<rmsd(2)
    rvec = fvec;
end

return;