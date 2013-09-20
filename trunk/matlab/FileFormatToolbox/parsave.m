function parsave(filename, varargin)
% parsave  Save workspace variables to file within a parfor loop.
%
% parsave(FILENAME, VAR1[, VAR2, ...])
%
%   Matlab's save() function cannot be used within a parfor loop of the
%   Parallel Computing Toolbox. This function bridges that gap.
%
%   FILENAME is a string with the path and filename to the file where the
%   variables will be saved to.
%
%   VAR1, VAR2, ... is a list of variables (note, the variables themselves,
%   not the name of the variables).
%
%   A comparison of save() vs. parsave(). Note the differences:
%
%     save('data.mat', 'tri', 'x')
%     parsave('data.mat', tri, x)
%
% This function is an extension of
% http://www.mathworks.co.uk/matlabcentral/fileexchange/30778-parsave
% to any number of variables.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
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
narginchk(2, Inf);
nargoutchk(0, 0);

% names of all input variables
for I = 2:nargin
    varname = genvarname(inputname(I));
    eval([varname ' = varargin{' num2str(I-1) '};'])
    
    % the first variable overwrites the previous file
    if (I == 2)
        save(filename, varname)
    else
        % other variables are added to the file
        save(filename, varname, '-append')
    end
    
end
