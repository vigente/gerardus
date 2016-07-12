% GRAB_SRC_CODE_DEPENDENCIES  Script to grab the Matlab functions required
% by the scripts used to run the experiments in this paper.
%
% Before running this script, we expect that:
%
%   * the script is in the paper's directory
%   "casero2015_3d_histology_diffusion_registration/src"
%
%   * the gerardus and private-gerardus toolboxes are in the path
%
% After running the script, directories will be created and Matlab
% functions copied over, following the rules:
%
%   * if the file is in ThirdPartyToolbox, we copy it and its licence file,
%   if available
%
%   * if the file is in a third party directory that corresponds to a whole
%   toolbox, we copy the whole toolbox, in order to e.g. preserve the
%   licence
%
%   * if the file is part of a gerardus toolbox, we simply copy that file

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2016 University of Oxford
% Version: 0.1.0
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Set Matlab's "Current Folder" to the directory where this script is
cd(fileparts(mfilename('fullpath')));

% list of top fscripts used in this paper
scripts = {'histology_blockface_mouse.m', ...
        'histology_noblockface_mouse.m', ...
        'histology_refinement_only_diffusion.m', ...
        'histology_refinement_without_diffusion.m', ...
        'manual_landmarks.m', ...
        'paper_figures.m', ...
        'toy_examples.m', ...
        'validate_reconstruction_hires.m', ...
        'validate_reconstruction.m', ...
        'validate_reconstruction_noblockface.m', ...
        'validate_reconstruction_onlydiff.m'};

% list all dependencies from all the scripts
fListTot = cell(0);
for file = scripts
    [fList, pList] = matlab.codetools.requiredFilesAndProducts(file);
    fListTot = [fListTot fList];
end

% remove duplicates
fListTot = unique(fListTot);

% remove top scripts from list of dependencies 
for I = 1:length(scripts)
    idx = strncmp(fListTot, pwd, length(pwd));
    fListTot(idx) = [];
end

% loop each dependency
for I = 1:length(fListTot)
    
    % find position of "/matlab/", which we are going to use to know where
    % the path to the gerardus or private-gerardus projects are
    idx = strfind(fListTot{I}, [filesep 'matlab' filesep]);
    
    if (isempty(idx))
        error(['String ''' filesep 'matlab' filesep ...
            ''' not found in I=' num2str(I) ': ' fListTot{I}])
    end
    
    % point to end of '/matlab/'
    idx = idx + length([filesep 'matlab' filesep]);
    
    % split the dependency into path, filename and extension
    [pth, file, ext] = fileparts(fListTot{I}(idx:end));

    % create directory to save the dependency if it doesn't exist
    if (isempty(dir(pth)))
        [ok, msg] = mkdir(pth);
        if (~ok)
            error(msg)
        end
    end
        
    % if the file is in ThirdPartyToolbox, we copy it and its licence file,
    % if available
    if strcmp(pth, 'ThirdPartyToolbox') % 'ThirdPartyToolbox' without subdirectory
        
        % copy the file
        [ok, msg] = copyfile(fListTot{I}, pth);
        if (~ok)
            error(msg)
        end
        
        % if the licence file is there, copy it too
        licfile = [fileparts(fListTot{I}), filesep, 'license-', file, '.txt'];
        if (exist(licfile, 'file') == 2)
            [ok, msg] = copyfile(licfile, pth);
            if (~ok)
                error(msg)
            end
        end
        
        licfile = [fileparts(fListTot{I}), filesep, 'licence-', file, '.txt'];
        if (exist(licfile, 'file') == 2)
            [ok, msg] = copyfile(licfile, pth);
            if (~ok)
                error(msg)
            end
        end
        
    % if the file is in a third party directory that corresponds to a whole
    % toolbox, we copy the whole toolbox, in order to e.g. preserve the
    % licence
    elseif ~isempty(strfind(pth, 'ThirdPartyToolbox'))
        
        % if we haven't copied this third party toolbox yet, we copy it
        if (length(dir('foo')) <= 2) % directory only contains "." and ".."
            [ok, msg] = copyfile(fileparts(fListTot{I}), pth);
            if (~ok)
                error(msg)
            end
        end
    
    % if the file is part of a gerardus toolbox, we simply copy that file
    else
        
        % copy the file
        [ok, msg] = copyfile(fListTot{I}, pth);
        if (~ok)
            error(msg)
        end
        
    end
    
end
