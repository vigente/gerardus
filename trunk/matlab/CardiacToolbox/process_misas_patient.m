function process_misas_patient(pathname)
% PROCESS_MISAS_PATIENT  Create MetaImage volumes from directory tree of
% IMA DICOM files, for project with human hearts provided by Kelvin Wong
%
% PROCESS_MISAS_PATIENT(PATHNAME)
%
%   PATHNAME is a string with the path to the parent directory containing
%   some directories with DICOM data. E.g.
%   PATHNAME='~/kelvin/MISAS_02_01_10'.

% Author: Ramon Casero <ramon.casero@comlab.ox.ac.uk>
% Copyright © 2009-2014 University of Oxford
% Version: 0.1.1
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% check arguments
narginchk(1, 1);
nargoutchk(0, 0);

% get list of all directories and files in the parent directory
dirlist = dir([pathname filesep '*']);

if isempty(dirlist)
    return
end

% for each data volume
for I = 1:length(dirlist)
    
    % skip files, we are interested in directories only
    if ~dirlist(I).isdir
        continue
    end
    
    % get a list of all files in that data volume, sorted correctly
    file = sortdirbynum(...
        [pathname filesep dirlist(I).name filesep '*.IMA'], 5, '.');
    
    if isempty(file)
        continue
    end
    
    % create a metaimage volume from the dicom files
    dcm2metaimage(...
        [pathname filesep dirlist(I).name filesep dirlist(I).name '.mha'], ...
        [], [], [], 'IMA', file);
    
    % move metaimage files to the parent directory
    movefile(...
        [pathname filesep dirlist(I).name filesep dirlist(I).name '.mha'], ...
        pathname)
    movefile(...
        [pathname filesep dirlist(I).name filesep dirlist(I).name '.raw'], ...
        pathname)
    
end

% if we have two NormalMri volumes, combine them with a max() function to
% remove banding artifacts, as suggested by Matt Robson

% check whether the files exist
fid1 = fopen([pathname filesep 'NormalMri1.mha'], 'r');
fid2 = fopen([pathname filesep 'NormalMri2.mha'], 'r');
if (fid1 ~= -1 && fid2 ~= -1)
    % close the files
    fclose(fid1);
    fclose(fid2);
    
    % load both files
    nrrd1 = scimat_load([pathname filesep 'NormalMri1.mha']);
    nrrd2 = scimat_load([pathname filesep 'NormalMri2.mha']);
    
    % to remove the banding artefacts, we take the maximum value for each voxel
    nrrd = nrrd1;
    nrrd.data = max(nrrd1.data, nrrd2.data);
    
    % save resulting file
    scimat_save([pathname filesep 'NormalMri.mat'], nrrd);
end

