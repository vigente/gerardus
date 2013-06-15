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

function statAnalysis(confs, objs1, objs2, method, cpath)

if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end
if (~exist([confs.OutDirectory '/Logs'],'dir'))
    mkdir([confs.OutDirectory '/Logs']);
end

switch method
    case 't_map'
        confs.OutputName = sprintf('%s_%s',confs.OutputNamePrefix,confs.Signal);
        diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' confs.OutputName '_t_stat.log']));
        % t-statistics
        performTstat(confs, objs1, objs2, cpath);
    case 'PCA'
        diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' confs.OutputName '_PCA_stat.log']));
        % PCA
        performPCA(confs, objs1);
end
diary('off');

return;
