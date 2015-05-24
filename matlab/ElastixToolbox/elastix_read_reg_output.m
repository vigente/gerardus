function [t, im, iterInfo] = elastix_read_reg_output(outdir)
% elastix_read_reg_output  Read output of registration computed with
% elastix.
%
% [T, IM, ITERINFO] = elastix_read_reg_output(OUTDIR)
%
%   OUTDIR is a string with the path to the output directory created by
%   elastix.
%
%   T is a struct with the contents of the parameter transform file
%   (OUTDIR/TransformParameters.0.txt). See elastix for details.
%
%   IM is the result image of the elastix registration (e.g.
%   OUTDIR/result.0.png, OUTDIR/result.0.jpg).
%
%   ITERINFO is a struct with the details of the elastix optimization
%   (OUTDIR/IterationInfo.0.R0.txt). See elastix for details.
%
% See also: elastix, blockface_find_frame_shifts.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
% Version: 0.2.2
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
nargoutchk(0, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TransformParameters.0.txt

t = elastix_read_file2param([outdir filesep 'TransformParameters.0.txt']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% result.0.png

% read result image
if (nargout > 1)
    
    % the result file can be in many different formats
    resultfile = dir([outdir filesep 'result.0.*']);
    
    if (isempty(resultfile))
        error(['No image result file: ' outdir filesep 'result.0.*'])
    end
    if (length(resultfile) > 1)
        error(['More than one result file: ' outdir filesep 'result.0.*'])
    end
    
    % read the image
    im = imread([outdir filesep resultfile.name]);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% IterationInfo.0.R0.txt

% read result image
if (nargout > 2)
    
    % the result file can be in many different formats
    iterfile = dir([outdir filesep 'IterationInfo.0.R0.txt']);
    
    if (isempty(iterfile))
        error(['No iteration info file: ' outdir filesep 'IterationInfo.0.R0.txt'])
    end
    if (length(resultfile) > 1)
        error(['More than one iteration info file: ' outdir filesep 'IterationInfo.0.R0.txt'])
    end
    
    % if the file is empty
    info = dir([outdir filesep 'IterationInfo.0.R0.txt']);
    if (info.bytes == 0)
        
        iterInfo.ItNr = [];
        iterInfo.Metric = [];
        iterInfo.StepSize = [];
        iterInfo.Gradient = [];
        iterInfo.Time = [];
        
    else
        
        % read the table
        table = dlmread([outdir filesep 'IterationInfo.0.R0.txt'], ...
            '\t', 1, 0);
        
        % create struct
        iterInfo.ItNr = table(:, 1);
        iterInfo.Metric = table(:, 2);
        iterInfo.StepSize = table(:, 3);
        iterInfo.Gradient = table(:, 4);
        iterInfo.Time = table(:, 5);
        
    end

end
