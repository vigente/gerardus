function get_gerardus_dependencies(files)
% GET_GERARDUS_DEPENDENCIES  Grab Gerardus files that a script depends on.
%
% This function processes one or more Matlab scripts to get a list of
% functions it depends on to run. It then creates a "gerardus" directory
% and makes a copy of those dependencies and other files that are needed to
% build MEX-files.
%
% Warning: This function is only a helper. In general it cannot create a
% subset of Gerardus that builds out of the box unless the script only
% depends on M-files. What it does automatically is that it grabs
%
%  * M-files
%  * MEX-files, their source files and CMakeLists.txt
%  * main gerardus directories (include, lib, programs, matlab, cpp)
%  * Files in the matlab directory, e.g. FindMatlab.cmake,
%    GerardusCommon.h, etc.
%
% Things that still need to be done manually:
%
%  * Copy C++ libraries that MEX-files depend on.
%  * Comment out CMakeLists.txt to remove unused MEX-files and libraries
%    from the build.
%  * Download and copy any external program not in the build, e.g. SCIP.
%
% GET_GERARDUS_DEPENDENCIES(FILES)
%
%   This function should be invoked from the same directory that contains
%   the scripts we want to process.
%
%   FILES is a cell array of strings. Each element is the name of a script
%   we want to process.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2015 University of Oxford
% Version: 0.1.0
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

% check arguments
narginchk(1, 1);
nargoutchk(0, 0);

% name of the local directories where we are going to put all the
% dependencies
GERLOCDIR = ['.' filesep 'gerardus'];
MATLOCDIR = [GERLOCDIR filesep 'matlab'];
INCLOCDIR = [GERLOCDIR filesep 'include'];
LIBLOCDIR = [GERLOCDIR filesep 'lib'];
PROGLOCDIR = [GERLOCDIR filesep 'programs'];
CPPLOCDIR = [GERLOCDIR filesep 'cpp'];
CPPSRCLOCDIR = [CPPLOCDIR filesep 'src'];
CPP3RDLOCDIR = [CPPSRCLOCDIR filesep 'third-party'];

%% copy files that are common to Gerardus

% path to Gerardus
% /home/jdoe/Software/gerardus/matlab
% /home/jdoe/Software/gerardus
matdir = fileparts(mfilename('fullpath'));
gerdir = fileparts(matdir);

% create common directories
mkdir_ifnotexist(GERLOCDIR);
mkdir_ifnotexist(INCLOCDIR);
mkdir_ifnotexist(LIBLOCDIR);
mkdir_ifnotexist(MATLOCDIR);
mkdir_ifnotexist(PROGLOCDIR);
mkdir_ifnotexist(CPP3RDLOCDIR);

% copy basic files
copyfile([gerdir filesep 'CMakeLists.txt'], GERLOCDIR);
copyfile([matdir filesep 'add_gerardus_paths.m'], MATLOCDIR);
copyfile([matdir filesep 'CMakeLists.txt'], MATLOCDIR);
copyfile([matdir filesep 'FindMatlab.cmake'], MATLOCDIR);
copyfile([matdir filesep 'LICENCE'], MATLOCDIR);
copyfile([matdir filesep 'GerardusCommon.h'], MATLOCDIR);
copyfile([matdir filesep 'GerardusCommon.hxx'], MATLOCDIR);
copyfile([matdir filesep 'MatlabExportFilter.h'], MATLOCDIR);
copyfile([matdir filesep 'MatlabExportFilter.hxx'], MATLOCDIR);
copyfile([matdir filesep 'MatlabImageHeader.h'], MATLOCDIR);
copyfile([matdir filesep 'MatlabImportFilter.h'], MATLOCDIR);
copyfile([matdir filesep 'MatlabImportFilter.hxx'], MATLOCDIR);
copyfile([matdir filesep 'MatlabMakeMacros.cmake'], MATLOCDIR);
copyfile([matdir filesep 'VectorWrapper.h'], MATLOCDIR);
copyfile([matdir filesep 'VectorWrapper.hxx'], MATLOCDIR);

%% copy all the M-files, MEX files, and CMakeLists.txt and source code 
%% files needed for the MEX-files

% parse a list of dependencies for first script
[fList, pList] ...
    = matlab.codetools.requiredFilesAndProducts(files);

% remove input scripts from the list of dependencies

% init variables
toolboxWithMexFile = cell(0);

% loop files
for I = 1:length(fList)
    
    % if the dependency is one of the input scripts, we skip it
    if any((cellfun(@(x)strcmp([pwd filesep x], fList{I}), files)))
        continue;
    end
    
    % name of the toolbox and file:
    % '/home/jdoe/Software/gerardus/matlab/ThirdPartyToolbox/iso2meshToolbox/mwpath.m'
    % 'ThirdPartyToolbox/iso2meshToolbox/mwpath.m'
    % {'ThirdPartyToolbox', 'iso2meshToolbox', 'mwpath.m'}
    toolboxname ...
        = regexprep(fList{I}, ['.*gerardus' filesep 'matlab' filesep], '');
    toolboxname = regexp(toolboxname, ['[^' filesep ']*'], 'match');
    
    % three possibilities:
    % * file is part of a third-party toolbox
    % * file is an M-file
    % * file is a MEX-file
    if strcmp(toolboxname{1}, 'ThirdPartyToolbox') % this file is part of a third-party toolbox
        
        % if this is a whole toolbox in the ThirdPartyToolbox directory
        if (length(toolboxname) == 3)
            
            % if we use a third party toolbox, we just copy the whole toolbox,
            % because we also want to include the licence; and in general it's not
            % a good idea to provide incomplete third-party software
            
            % create directory for toolbox
            itExisted = mkdir_ifnotexist([MATLOCDIR filesep toolboxname{1} filesep toolboxname{2}]);
            
            % check that we have not already copied the third party toolbox
            if (~itExisted)
                
                % copy third party toolbox content
                ok = copyfile(fileparts(fList{I}), ...
                    [MATLOCDIR filesep 'ThirdPartyToolbox' filesep toolboxname{2}]);
                if (~ok)
                    error(['Cannot copy files to directory: '...
                        MATLOCDIR filesep 'ThirdPartyToolbox' filesep toolboxname{2}])
                end
                
            end
            
        elseif (length(toolboxname) == 2) % this is only a function in the ThirdPartyToolbox directory
            
            % create the ThirdPartyToolbox if it doesn't exist
            mkdir_ifnotexist([MATLOCDIR filesep 'ThirdPartyToolbox']);
            
            % copy the file
            copyfile([matdir filesep 'ThirdPartyToolbox' filesep toolboxname{2}], ...
                [MATLOCDIR filesep 'ThirdPartyToolbox'])
            
            % copy the licence file, if available
            [~, filename, ext] = fileparts([matdir filesep 'ThirdPartyToolbox' filesep toolboxname{2}]);
            licfile = [matdir filesep 'ThirdPartyToolbox' filesep 'licence-' filename '.txt'];
            if (exist(licfile, 'file'))
                copyfile(licfile, [MATLOCDIR filesep 'ThirdPartyToolbox']);
            end
            
        end
            
    elseif (exist(fList{I}, 'file') == 2) % not 3rd party, this is an M-file
        
        % create directory for toolbox if necessary
        mkdir_ifnotexist([MATLOCDIR filesep toolboxname{1}]);
        
        % copy the M-file
        ok = copyfile(fList{I}, ...
            [MATLOCDIR filesep toolboxname{1} filesep toolboxname{2}]);
        if (~ok)
            error(['Cannot copy files to directory: '...
                MATLOCDIR filesep 'ThirdPartyToolbox' filesep toolboxname{2}])
        end
        
    elseif (exist(fList{I}, 'file') == 3) % not 3rd party, this is a MEX file
        
        % get path and file name
        [pth, mexname, ext] = fileparts(fList{I});
        
        % full path to CMakeLists.txt
        cmakelistspth = [pth filesep 'CMakeLists.txt'];
        
        % check that a CMakeLists.txt exists for this toolbox
        if (~exist(cmakelistspth, 'file'))
            error(['No CMakeLists.txt found for function: ' ...
                fList{I}])
        end
        
        % copy CMakeLists.txt
        ok = copyfile([pth filesep 'CMakeLists.txt'], [MATLOCDIR filesep toolboxname{1}]);
        if (~ok)
            error(['Cannot copy file: ' pth filesep 'CMakeLists.txt'])
        end
        
        % get the names of the source files that we need to build the mex file
        src = grepCodeFile(cmakelistspth, mexname);
        
        % create directory for toolbox if necessary
        mkdir_ifnotexist([MATLOCDIR filesep toolboxname{1}]);
        
        % copy the source files
        for J = 1:length(src)
            ok = copyfile([pth filesep src{J}], [MATLOCDIR filesep toolboxname{1}]);
            if (~ok)
                error(['Cannot copy file: ' src{J}])
            end
        end
        
        % keep track of the toolboxes that contain MEX files (without
        % repetition)
        toolboxWithMexFile(end+1) = {toolboxname{1}};
        toolboxWithMexFile = unique(toolboxWithMexFile);
        
    else
        
        error('Function is neither an M- nor a MEX-file')
        
    end

end


end

% create directory if it doesn't exist 
function itExisted = mkdir_ifnotexist(pth)

% check whether the directory already exists
itExisted = exist(pth, 'dir');

% if toolbox directory doesn't exist
if (~itExisted)
    
    % create toolbox directory
    ok = mkdir(pth);
    if (~ok)
        error(['Cannot create directory: ' pth])
    end
    
end

end

% returns a cell array with the source code files necessary to compile
% target s in CMakeLists.txt (file is the path and full file name)
function sout = grepCodeFile(file, s)

% read file into a string
sfile = fileread(file);

% grep the bit where the mex target is defined
% add_mex_file(cgal_check_self_intersect
%  CgalCheckSelfIntersect.cpp
sout = regexpi(sfile, ['add_mex_file\(' s '[^\)]*'], 'match');
sout = sout{1};

% split into words
% 'add_mex_file(cgal_check_self_intersect'    'CgalCheckSelfIntersect.cpp'    ''
sout = strsplit(sout);

% remove first and empty words
sout = sout(2:end);
sout(cellfun(@isempty, sout)) = [];

end

