function [ir, pr, jc] = sparse_breakdown(s)
% SPARSE_BREAKDOWN  Extract internal arrays from sparse matrix
%
% [IR, PR, JC] = SPARSE_BREAKDOWN(S)
%
%   S is a sparse matrix.
%
%   IR is the vector obtained with the C++ function mxGetIr().
%   http://www.mathworks.com/access/helpdesk/help/techdoc/apiref/mxgetir.html
%
%   PR is the vector obtained with the C++ function mxGetPr().
%   http://www.mathworks.com/access/helpdesk/help/techdoc/apiref/mxgetpr.html
%
%   JC is the vector obtained with the C++ function mxGetJc().
%   http://www.mathworks.com/access/helpdesk/help/techdoc/apiref/mxgetjc.html

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010 University of Oxford
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
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% It is expected that sparse_breakdown.cpp has been compiled with mex:
%
% 32 bit
%   >> mex sparse_breakdown.cpp
%
% 64 bit
%   >> mex -largeArrayDims sparse_breakdown.cpp
error('Compiled MEX function has not been found')
