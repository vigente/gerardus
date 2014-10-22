function [y] = JoFunction(x, p)
% JoFunction  Power of a scalar.
%
% Y = myfun(X, P)
% 
%   X is a scalar.
%
%   Y is X^P.

% Author: Jo Bates <jobates81@gmail.com>
% Copyright © 2014 University of Oxford
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

% check arguments
narginchk(2, 2);
nargoutchk(0, 1);

% defaults
if (nargin < 2 || isempty(p))
    p = 1;
end
if (~isscalar(x))
  error('X must be a scalar.')
end
y = x.^p;
end

