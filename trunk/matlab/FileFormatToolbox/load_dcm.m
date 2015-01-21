function [dcm, dcminfo] = load_dcm(files, pth)
% LOAD_DCM  Load array of DICOM filenames as 3D+t image and metainformation
% matrix.
%
% [DCM, DCMINFO] = LOAD_DCM(FILES, PTH)
%
%   FILES, PTH are the output variables of function load_multi_dcm_dir().
%
%   FILES is an array of structs with information for each DICOM file in
%   the set. FILES(I) corresponds to the I-th slice. Each slice contains
%   files for each temporal frame.
%
%   PTH is a string with the path to the DICOMDIR file.
%
%   DCM is a 4-D array with the DICOM images. The format is 
%   DCM(row, column, slice, time frame).
%
%   DCMINFO is a matrix of structs. DCMINFO(slice, time frame) contains the
%   metainformation for the corresponding slice and time frame in the data
%   set.
%
% See also: load_multi_dcm_dir.

% Authors: Ramon Casero <rcasero@gmail.com> and Benjamin Villard
% <b.016434@gmail.com>
% Copyright © 2014-2015 University of Oxford
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
% http://code.google.com/p/gerardus/
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
narginchk(1, 2);
nargoutchk(0, 2);

% defaults
if (nargin < 2)
    pth = '.';
end

% loop list of files (slices)
for S = 1:length(files)
    
    % loop frames within each file
    for T = 1:length(files(S).ImageNames)
        
        % load image metainformation. We switch the order of the indices
        % because I corresponds to temporal frames and J corresponds to
        % physical slices. In the image output, we are also going to switch the
        % time/slice indices the same way
        if ((T == 1) && (S == 1))
            % when we read the metainformation from the first frame in the
            % first file, we allocate memory for all the other metainformation
            % structs
            dcminfo = dicominfo([pth filesep files(S).ImageNames{T}]);
            dcminfo = repmat(dcminfo, length(files), length(files(S).ImageNames));
        else
            dcminfo(S, T) = dicominfo([pth filesep files(S).ImageNames{T}]);
        end
        
        % load images
        if ((T == 1) && (S == 1))
            % when we read the first frame from the first image, we allocate
            % memory for all the other frames
            dcm = dicomread(dcminfo(1, 1));
            dcm = repmat(dcm, 1, 1, length(files), length(files(S).ImageNames));
            
        else
            
            dcm(:, :, S, T) = dicomread(dcminfo(S, T));
            
        end
        
    end
    
end

