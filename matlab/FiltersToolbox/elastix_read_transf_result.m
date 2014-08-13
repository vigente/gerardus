function [param, im, iterInfo] = elastix_read_transf_result(outdir)
% elastix_read_transf_result  Read the result image and transform
% parameters of a registration computed with elastix
%
% [PARAM, IM, ITERINFO] = elastix_read_transf_result(OUTDIR)
%
%   OUTDIR is a string with the path to the output directory created by
%   elastix.
%
%   PARAM is a struct with the contents of the parameter transform file
%   (OUTDIR/TransformParameters.0.txt).
%
%   IM is the result image of the elastix registration (e.g.
%   OUTDIR/result.0.png, OUTDIR/result.0.jpg).
%
%   ITERINFO is a struct with the details of the elastix optimization
%   (OUTDIR/IterationInfo.0.R0.txt).

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
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
nargoutchk(0, 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% TransformParameters.0.txt

% init output struct
param = struct([]);

% open transform parameters text file to read
fid = fopen([outdir filesep 'TransformParameters.0.txt'], 'r');
if (fid == -1)
    error(['Cannot open ' outdir filesep 'TransformParameters.0.txt'])
end

% read text file line by line
tline = fgetl(fid);
while ischar(tline)
    
    % duplicate line to process it
    tline2 = tline;
    
    % skip lines that have no parameters
    if (isempty(tline2) || (tline2(1) ~= '('))
        
        % read next line
        tline = fgetl(fid);
        
        continue;
    end
    
    % sanity check that parameter line ends in ')'
    if (tline2(end) ~= ')')
        error(['Line does''t end in '')'': ' tline2])
    end
    
    % remove '(' and ')' from the line
    tline2 = tline2(2:end-1);
    
    % position of first whitespace
    idx = strfind(tline2, ' ');
    
    % get parameter field label
    if (isempty(idx))
        label = tline2;
    else
        label = tline2(1:idx-1);
    end
    
    % get parameter values
    if (isempty(idx))
        val = [];
    else
        val = tline2(idx+1:end);
    end
    
    % if the value can be converted to numeric format, we do. Otherwise,
    % it's a string, and we remove the double quotes ""
    val2 = str2num(val);
    if (~isempty(val2))
        val = val2;
    else
        if (val(1) == '"')
            val = val(2:end);
        end
        if (val(end) == '"')
            val = val(1:end-1);
        end
    end
    
    % add field to output struct
    if (isempty(param))
        param = struct(label, val);
    else
        param.(label) = val;
    end
    
    % read next line
    tline = fgetl(fid);
    
end

% close parameters file
fclose(fid);

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
    
    % read the table
    table = dlmread([outdir filesep 'IterationInfo.0.R0.txt'], '\t', 1, 0);
    
    % create struct
    iterInfo.ItNr = table(:, 1);
    iterInfo.Metric = table(:, 2);
    iterInfo.stepSize = table(:, 3);
    iterInfo.Gradient = table(:, 4);
    iterInfo.Time = table(:, 5);
    
end

