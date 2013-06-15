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

function outNames = SpharmMatExpansion(confs, objs, method)

numSbj = length(objs);
if ~exist(confs.OutDirectory,'dir')
    mkdir(confs.OutDirectory);
end
if (~exist([confs.OutDirectory '/Logs'],'dir'))
    mkdir([confs.OutDirectory '/Logs']);
end
outNames = {};

h = waitbar(0,'Please wait...');
for i = 1:numSbj
    file = objs{i};
    [path, name, ext] = fileparts(file);
    
    switch method
        case 'ExpLSF'
            diary(fullfile([confs.OutDirectory '/Logs'],[date, '_', num2str(round(cputime)), '_' name '_ExpLSF.log']));
            % Call LSF method
            maxDeg = confs.MaxSPHARMDegree;
            [fvec, deg, Z, outNames{end+1}]=create_SPHARM_des_LSF([],[],[],maxDeg,file,deblank(char(confs.OutDirectory)));
            clear('fvec','deg','Z');
    end
    diary('off');
    waitbar(i/numSbj)
end
close(h);

return