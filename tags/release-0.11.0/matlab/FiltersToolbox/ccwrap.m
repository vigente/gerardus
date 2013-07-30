function varargout = ccwrap(f, varargin)
% CCWRAP  Wrapper that enables CTRL+C interruption of MEX files
%
% Matlab code (.m) can be interrupted by the user pressing the keys CTRL+C.
% This is not the case with MEX functions (e.g. C++ code compiled to run as
% a MEX function).
%
% To overcome this limitation, this function can be used as a wrapper that
% enables CTRL+C interruption of MEX functions.
%
% Note that this wrapper prevents memory leaks. Even if your C++ code
% allocates memory using functions other than Matlab's
% mxCreateDoubleMatrix, mxCreateNumericArray, etc, it will run your MEX
% function on  a worker from the Parallel Toolbox. When the computations
% finish (or CTRL+C is pressed), the worker is closed, and all memory is
% released.
%
% VARARGOUT = ccwrap(@F, VARARGIN)
%
%   @F is a function handle to the MEX function.
%
%   VARARGIN is the list of input arguments to function F.
%
%   VARARGOUT is the list of output arguments from function F.
%
% Example:
%
%   % this MEX function cannot be interrupted pressing CTRL+C
%   im2 = itk_imfilter('bwdilate', 3, 1);
%
%   % this MEX function now can be interrupted pressing CTRL+C
%   im2 = ccwrap(@itk_imfilter, 'bwdilate', 3, 1);
%
% See also: C++ function ctrlcCheckPoint() in 
% http://code.google.com/p/gerardus/source/browse/trunk/matlab/GerardusCommon.hpp


% Authors: Ramon Casero <rcasero@gmail.com>, Friedrich Hempel (MathWorks).
% Copyright © 2012 University of Oxford
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
narginchk(1, Inf);
nargoutchk(0, Inf);

% f must be a handle to a function
if ~isa(f, 'function_handle')
    error('F must be a function handle');
end

% set up a job that will run input function f as a task
%   note: onCleanup() needs an output argument, even if we don't use it.
%   Otherwise, the cleanup function will execute immediately, instead of
%   waiting for an exit
c = parcluster;
job = createJob(c);
foo = onCleanup(@() cleanup(job, f));
 
% launch a task to run input function f
createTask(job, f, nargout, varargin);
submit(job);
wait(job);

% get results
varargout = fetchOutputs(job);
 
end

% cleanup(): function called when ccwrap() exits. This can be due to two
% reasons: 1) the task finished successfully or 2) the user interrupted
% ccwrap() with CTRL+C
function cleanup(job, f)

switch job.State
    case 'running'
        job.cancel;
        aux = functions(f);
        [pathstr, name] = fileparts(aux.file);
        fprintf(2, ['Operation terminated by user during ' ...
            '<a href="matlab:helpUtils.errorDocCallback(''%s'',''%s%s%s.m'',0)">' ...
            '%s</a>\n'], aux.function, pathstr, filesep, name, aux.function);
    case 'finished'
        % DEBUG:
%         disp('my_cleanup called because job finished');
end

end
