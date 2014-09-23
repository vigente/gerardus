function [tinv, iminv] = elastix_invtransf(t, opts)
% elastix_invtransf  Compute inverse of transform for elastix.
%
% [TINV, IMINV] = elastix_invtransf 

% Author: Ramon Casero <rcasero@gmail.com>
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
narginchk(1, 2);
nargoutchk(0, 2);

% defaults
if (nargin < 2)
    opts = [];
end

switch (t.Transform)
    
    case 'BSplineTransform'
        
        if (~isfield(opts, 'regParam'))
            error('For BSplineTransform, opts.regParam must be provided')
        end

        % set transform that we will use to invert the input transform
        opts.regParam.HowToCombineTransforms = 'Compose';
        opts.regParam.Metric = 'DisplacementMagnitudePenalty';
        
        % set the transform we want to invert as the initial transform
        el_opts.t0 = t;
        el_opts.verbose = 1;
        
        % invert the transform, as described in the elastix manual and in Metz et
        % al. [2011]
        [tinv, iminv] = elastix(opts.regParam, opts.im, opts.im, el_opts);

    otherwise
        
        % TODO: affine transforms can be inverted easily by centering and
        % then inverting the affine matrix, but I have no time now
        
        error('Transform not implemented')
        
end
