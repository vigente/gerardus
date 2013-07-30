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
% run shrec for one iteration
%

function [fvec,rmsd] = shrec_kth_itr(fvec,atlas,max_d,k);
global alpha;
global beta;
global gamma;
global gran;
global res;
% global meshsize;
global dg;

rmsd(1) = SPHARM_rmsd(fvec, atlas);

if rmsd(1)==0
    rmsd = [0 0 0]; return;
end

if k==1
    % in the first iteration, fvec and atlas may not be aligned in either
    % object or parameter space
    
    % option 1: assume parameter space not aligned
    [fvec_a, Robj] = match_icp(fvec, atlas, max_d); % object space
    [fvec_a, Aprm] = match_param_hie(fvec_a,atlas,alpha,beta,gamma,gran,res,max_d); % parameter space
    [fvec_a, Robj] = match_cps(fvec_a, atlas, max_d);
    rmsd_a = SPHARM_rmsd(fvec_a, atlas);
    
    % option 2: assume parameter space aligned
    [fvec_b, Robj] = match_cps(fvec, atlas, max_d);
    rmsd_b = SPHARM_rmsd(fvec_a, atlas);
    
    % pick the better one
    if rmsd_a > rmsd_b
        fvec = fvec_b; rmsd(2) = rmsd_b;
    else
        fvec = fvec_a; rmsd(2) = rmsd_a;
    end    
        
else
    % assume aligned in parameter space, refine registration in the object space
    [fvec, Robj] = match_cps(fvec, atlas, max_d);
    rmsd(2) = SPHARM_rmsd(fvec, atlas);
end

% assume aligned in object space, register the parameterization
[fvec, Aprm] = match_param_hie(fvec,atlas,alpha,beta,gamma,gran,res,max_d);
rmsd(3) = SPHARM_rmsd(fvec, atlas);

disp(sprintf('Base Res %d, Step %d, Depth %d, Top %d, Gamma %d: RMSD %f => %f => %f',res,rmsd));
    
return;
