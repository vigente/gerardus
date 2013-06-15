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
% distance between two SPHARM surfaces
%

function rmsd = SPHARM_rmsd(fvec1, fvec2)

dist = fvec1-fvec2;
rmsd = norm(dist(:))/sqrt(4*pi);

return;

%
% distance between two SPHARM surfaces using sampling
%

function rmsd = PDM_rmsd(fvec1, fvec2, res)

vn = size(fvec1,1);
max_d = sqrt(vn)-1;

% res = -3
[vs, fs] = sample_Param_Space(res);
Z = calculate_SPHARM_basis(vs, max_d);

vs1 = real(Z*fvec1);
vs2 = real(Z*fvec2);

rmsd = (vs1-vs2).^2;
rmsd = sqrt(sum(rmsd(:))/size(vs,1));

return;