function [dcm, dcminfo] = load_dcm(files, pth)
%% LOAD_DCM  Creates two 4-D arrays from DICOM series Data files.
%
% LOAD_DCM reads the data paths acquired via Christopher Kelly's LOAD_MULTI_DCM_DIR function. 
% The data is arranged in a 2D array containing the dicominfo (metainfo) for each individual
% slice across time, and a 4D array containing the slices (.dcm) with the
% rows, columns, slice and time.
%
% "S" represents the Slices, and "T" represents the Frames, which are
% temporal data. 
% 
% The dcm output has for format: dcm(r,c,S,T).
% The dcminfo output has for format: dcminfo(S,T)
% Where "r" and "c" are the rows and columns respectively of the slice "S"
% belonging to frame "T". 
%
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

