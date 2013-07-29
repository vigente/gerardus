% ADD_GERARDUS_PATHS  Add all Gerardus toolboxes to the Matlab path
%
%   add_gerardus_paths
%
%   This script takes no input arguments and produces no output arguments.
%   It will check the current Matlab path, and add to it those Gerardus
%   toolboxes that are not already in the path.
%
%   It then asks the user whether the updated path should be saved to
%   gerardus/matlab/pathdef.m.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
% Version: 0.2.0
% $Rev$
% $Date$
%
% University of Oxford means the Chancellor, Masters and Scholars of
% the University of Oxford, having an administrative office at
% Wellington Square, Oxford OX1 2JD, UK. 
%
% This file is part of Gerardus.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. The offer of this
% program under the terms of the License is subject to the License
% being interpreted in accordance with English Law and subject to any
% action against the University of Oxford being under the jurisdiction
% of the English Courts.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

% path to the gerardus/matlab directory (or equivalently, path to this file
% add_gerardus_paths.m)
matdirpath = mfilename('fullpath');
matdirpath = fileparts(matdirpath);

% list of toolboxes to add to the path
toolboxes = { ...
    'CardiacToolbox', ...
    'CgalToolbox', ...
    'FileFormatToolbox', ...
    'FiltersToolbox', ...
    'ItkToolbox', ...
    'ManifoldToolbox', ...
    'PointsToolbox', ...
    'ThirdPartyToolbox', ...
    ['ThirdPartyToolbox' filesep 'FastMarchingToolbox'], ...
    ['ThirdPartyToolbox' filesep 'FastMarchingToolbox/toolbox'], ...
    ['ThirdPartyToolbox' filesep 'iso2meshToolbox'], ...
    ['ThirdPartyToolbox' filesep 'iso2meshToolbox/bin'], ...
    ['ThirdPartyToolbox' filesep 'niftiToolbox'], ...
    ['ThirdPartyToolbox' filesep 'SphereUniformSamplingToolbox'], ...
    ['ThirdPartyToolbox' filesep 'SphericalTrigToolbox'], ...
    ['ThirdPartyToolbox' filesep 'SpharmToolbox/code'], ...
    ['ThirdPartyToolbox' filesep 'SpharmToolbox/code/C_sources'], ...
    ['ThirdPartyToolbox' filesep 'sphsplineToolbox'], ...
    ['ThirdPartyToolbox' filesep 'TiffreadToolbox'] ...
    };

% current Matlab paths
pth = path;

for I = 1:length(toolboxes)
    
    % if the toolbox is not in the path already, add it
    if isempty(strfind(pth, [matdirpath filesep toolboxes{I}]))
        disp(['Adding ' toolboxes{I} ' to the path'])
        
        addpath([matdirpath filesep toolboxes{I}]);
    else
        disp(['Skipping ' toolboxes{I} '. Already in the path'])
    end
    
end

% ask user whether to save the new path
str = input('Save current path to gerardus/matlab (y/N)? ', 's');
if (str == 'y')
    savepath([matdirpath filesep 'pathdef.m'])
end

%% System paths to DLLs for Windows
if ~isempty(strfind(getenv('OS'), 'Windows'))
    
    % full paths to directories with DLLs
    pathsToDlls = {cd(cd('..\lib')), cd(cd('..\lib\bin')), cd(cd('..\cpp\src\third-party\CGAL-4.2\auxiliary\gmp\lib'))};
    
    % get system paths
    systemPath = getenv('PATH');
    
    % add paths to directories with DLLs unless they are already in the
    % system path
    for I = 1:length(pathsToDlls)
        if isempty(strfind(lower(systemPath), lower(pathsToDlls{I})))
            disp(['Adding ' pathsToDlls{I} ' to system path'])
            systemPath = [pathsToDlls{I} ';' systemPath];
        end
    end
    setenv('PATH', systemPath);
end
