function [tout, info] = transfdiff(opt, tp, tm)
% TRANSFDIFF Transform diffusion algorithm for sequence of images.
%
% This function is part of a larger algorithm that solves the following
% registration problem:
%
% We have a sequence of images I=1,2,...,N. We want to register each image
% to its adjacent neighbours to form a volume. But instead of registering
% I=2 to I=1, and then I=3 to the result, and then I=4 to the result, etc.,
% it can be shown that all images can be aligned in parallel using an
% iterative process that we call "transform diffusion registration".
%
% TRANSFDIFF solves the "transform diffusion" step in the algorithm:
%
% Let's assume that somebody has already registered each image I to its
% neighbours I-1 and I+1. TRANSFDIFF takes that set of registration
% parameters and "diffuses" them to produce a transform for each image.
% Applying these transforms to the images produces the globally aligned
% sequence.
%
%
% TRANSFDIFF applies a diffusion process to those neighbour transforms. The
% output is the N transforms that applied to the N images best align them.
%
% [TOUT, INFO] = TRANSFDIFF(OPT, ...)
%
%   OPT is a string with the name of the transform to apply, or a struct
%   with the transform name and algorithm parameters. The syntax of the
%   function depends on the transform (see below).
%
%   TOUT is an array with N output transforms for the N slices.
%
%   INFO is a struct with information about the algorithm internals:
%
%     'NumIter': Number of diffusion iterations until stop condition.
%
% -------------------------------------------------------------------------
%
% * OPT='TranslationTransform'
% * OPT.Transform='TranslationTransform'
%
% [TOUT, INFO] = TRANSFDIFF(OPT, TP)
%
%   TP is a matrix with N-1 rows. TP(I, :) is the translation from slice I
%   to slice I+1.
%
%   This case assumes that the transforms are symmetric, i.e. the
%   translation from slice I+1 to I is -TP(I, :).
%
% [TOUT, INFO] = TRANSFDIFF(OPT, TP, TM)
%
%   TM is similar to TP, but TM(I, :) is the translation from slice I+1 to
%   I. This syntax is used when the transform is non-symmetric (e.g. when
%   TP, TM are B-spline coefficients).
%
% -------------------------------------------------------------------------
%
% * OPT='AffineTransform'
% * OPT.Transform='AffineTransform'
%
% * OPT='EulerTransform'
% * OPT.Transform='EulerTransform'
%
% [TOUT, INFO] = TRANSFDIFF(OPT, TP)
%
%   TP is an array of affine transforms in homogeneous coordinates. TP(I,:)
%   is the affine transform from slice I to slice I+1. Each transform has
%   the format
%
%   TP(:, :, I) = [A 0]
%                 [d 1]
%
%   where A is a matrix and d is a translation vector, and the affine
%   transform of a row vector X is defined as Y=X*TP(:,:,I).
%
%   This case assumes that the transforms are symmetric, i.e. the
%   translation from slice I+1 to I is inv(TP(I, :)).
%
% [TOUT, INFO] = TRANSFDIFF(OPT, TP, TM)
%
%   As in the translation case, TM is similar to TP, but TM(I, :) is the
%   translation from slice I+1 to I. This syntax is used when the transform
%   is non-symmetric
%
% See also: transdiffreg.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2015 University of Oxford
% Version: 0.2.4
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
narginchk(2, 3);
nargoutchk(0, 2);

% defaults
if (isempty(opt))
    error('OPT cannot be empty');
elseif (ischar(opt))
    aux = opt;
    clear opt
    opt.Transform = aux;
end
if (~isfield(opt, 'Epsilon'))
    opt.Epsilon = 1e-3;
end
if (~isfield(opt, 'MaxIter'))
    opt.MaxIter = Inf;
end
if (~isfield(opt, 'Alpha'))
    opt.Alpha = 0.49;
end

% convert to array format if transforms are provided in elastix format
[tp, tp0] = elastix2mat(tp);
if (nargin > 2)
    [tm, tm0] = elastix2mat(tm);
end

% preprocessing of inputs
switch (opt.Transform)
    
    case {'TranslationTransform', 'BSplineTransform'}
        
        % compute "tm" from "tp"
        if (nargin < 3 || isempty(tm))
            tm = -tp;
        end
        
        % add dummy transforms at the beginning of "tm" and at the end of
        % "tp". This makes it easier to operate, because then:
        % * tp(I) and tm(I) are the neighbours of I for intermediate slices
        % * extreme slices can be dealt with in the same way as
        %   intermediate slices
        tp = tp([1:end end], :);
        tm = tm([1 1:end], :);
        tp(end, :) = tm(end, :);
        tm(1, :) = tp(1, :);
        
    case {'EulerTransform', 'AffineTransform'}
        
        % compute "tm" from "tp"
        if (nargin < 3 || isempty(tm))
            for I = 1:size(tp, 3)
                tm(:, :, I) = inv(tp(:, :, I));
            end
        end
        
        % add dummy transforms at the beginning of "tm" and at the end of
        % "tp". This makes it easier to operate (see note above)
        tp = tp(:, :, [1:end end]);
        tm = tm(:, :, [1 1:end]);
        tp(:, :, end) = tm(:, :, end);
        tm(:, :, 1) = tp(:, :, 1);
        
    otherwise
        
        error(['Transform ' opt.Transform ' not implemented'])
        
end

% check inputs
if any(size(tp) ~= size(tm))
    error('TP and TM must have the same size')
end

% if no transforms provided, then we don't need to run the algorithm
if (isempty(tp))
    tout = [];
    return
end

% apply registration diffusion
switch (opt.Transform)
    
    case {'TranslationTransform', 'BSplineTransform'}
        
        [tout, info] = translation_diffusion(opt, tp, tm);

    case {'EulerTransform', 'AffineTransform'}
        
        [tout, info] = affine_diffusion(opt, tp, tm);
        
    otherwise
        
        error('Transform not implemented')
        
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOCAL FUNCTIONS

%% convert transform from elastix format to array, if necessary
function [tout, t0] = elastix2mat(t)

if (isstruct(t))
    
    % save copy of input so that we can convert the output back to elastix
    % format
    t0 = t;
    
    switch (t(1).Transform)
        
        case {'TranslationTransform', 'BSplineTransform'}
            
            tout = cat(1, t(:).TransformParameters);
            
        % TODO: AffineTransform
            
        otherwise
            
            error('Conversion from this elastix transform to Matlab array not implemented')
    end
    
else % already in array format
    
    % no need to save a copy of the input
    t0 = [];
    
    % output is just the input
    tout = t;
    
end

end

%% translation_diffusion: Apply diffusion to translation transforms
function [tout, info] = translation_diffusion(opt, tp, tm)

% init accumulated transform to apply to each slice
tout = zeros(size(tp));

I = 0;
while (1)

    % iteration number
    I = I + 1;

    % transform to apply to each slice in this iteration (transform update)
    t = opt.Alpha * (tp + tm);
    
    % compose with the previous transforms to obtain a total transform
    tout = tout + t;
    
    % update neighbour transforms ...
    % ... of internal images
    tp(1:end-1, :) = t(2:end, :) + tp(1:end-1, :) - t(1:end-1, :);
    tm(2:end, :) = t(1:end-1, :) + tm(2:end, :) - t(2:end, :);
    
    %... of extreme images
    tp(end, :) = tm(end, :);
    tm(1, :) = tp(1, :);
    
    % stopping criteria...
    % ... norm of the transform update
    tnorm = sqrt(sum(t.^2, 2));
    info.MaxTNorm(I) = max(tnorm);
    info.MeanTNorm(I) = mean(tnorm);
    if (all(tnorm <= opt.Epsilon))
        break
    end
    
    % ... maximum number of iterations
    if (I == opt.MaxIter)
        break
    end
    
end

info.NumIter = I;

end

%% affine_diffusion: Apply diffusion to affine transforms
function [tout, info] = affine_diffusion(opt, tp, tm)

% init accumulated transform to apply to each slice
tout = eye(size(tp(:, :, 1)));
tout = repmat(tout, 1, 1, size(tp, 3));

% allocate memory for transform update
t = zeros(size(tp));

I = 0;
while (1)

    % iteration number
    I = I + 1;
    
    % compute the transform to apply to each slice in this iteration
    for J = 1:size(tp, 3)
        
        % transform to apply to each slice in this iteration (transform
        % update)
        t(:, :, J) = real(expm(opt.Alpha ...
            * (logm(tp(:, :, J)) + logm(tm(:, :, J)))));
        
        % compose with the previous transforms to obtain a total transform
        tout(:, :, J) = tout(:, :, J) * t(:, :, J);
        
    end
    
    % update neighbour transforms of internal images
    for J = 2:size(tm, 3)
        tm(:, :, J) = t(:, :, J) \ tm(:, :, J) * t(:, :, J-1);
    end
    for J = 1:size(tp, 3)-1
        tp(:, :, J) = t(:, :, J) \ tp(:, :, J) * t(:, :, J+1);
    end
        
    % update neighbour transforms of extreme images
    tp(:, :, end) = tm(:, :, end);
    tm(:, :, 1) = tp(:, :, 1);
    
    % stopping criteria...
    % ... norm of the transform update different from identity
    tnorm = zeros(1, size(t, 3));
    for J = 1:size(t, 3)
        tnorm(J) = norm(t(:, :, J) - eye(size(t(:, :, 1))), 2);
    end
    info.MaxTNorm(I) = max(tnorm);
    info.MeanTNorm(I) = mean(tnorm);
    if (all(tnorm <= opt.Epsilon))
        break
    end
    
    % ... maximum number of iterations
    if (I == opt.MaxIter)
        break
    end
    
end

info.NumIter = I;

end
