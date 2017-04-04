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
%   with the transform name and algorithm parameters. Options common to all
%   transformations:
%
%     'MaxIter': (def 50) Stopping criterion. The algorithm stops after
%                MaxIter diffusion iterations.
%
%     'Alpha':   (def 0.45) Diffusion coefficient. The smaller the value of
%                alpha, more iterations are needed for the same result. But
%                alpha >0.45 may not smooth high frequencies very well, and
%                values >=0.5 cause oscillations.
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
%   Specific OPT options:
%
%     'Epsilon': (def 0.0) Stopping criterion. It stops when all
%                translation components are small:
%                for all t, abs(t)<=Epsilon.
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
%   transformation from slice I+1 to I is inv(TP(I, :)).
%
%   Specific OPT options:
%
%     'Epsilon': (def 0.0) Stopping criterion. It stops when the norm of
%                the difference between the transformation and the identity
%                matrix is small:
%                for all t, norm(t-eye(3))<=Epsilon.
%
% -------------------------------------------------------------------------
%
% * OPT.Transform='BsplineTransform'
%
% [TOUT, INFO] = TRANSFDIFF(OPT, TP, TM)
%
%   TP is a matrix with N-1 rows. Each row contains the B-spline
%   coefficients with the transformation from slice I to slice I+1. The
%   format of the B-spline coefficient vector is 
%   [p0x, p1x, ..., pNx, p0y, p1y, ..., pNy].
%
%   TM is similar to TP, but TM(I, :) is the transformation from slice I+1
%   to I.
%
%   Specific OPT options:
%
%     'Epsilon': (def 0.0) Stopping criterion. It stops when all
%                coefficient displacement components are small:
%                for all t, abs(t)<=Epsilon.
%
%     'tbsp': (def []) For Choi and Lee (2000) injectivity criteria. These
%             are sufficient criteria that guarantee the B-spline doesn't
%             fold over. Transform struct in elastix format. The
%             tbsp.TransformParameters are ignored, only the grid
%             information is used. If not provided or empty, the conditions
%             are not used as a stopping criterion, and the resulting
%             B-splines may fold over.
%
%
% -------------------------------------------------------------------------
%
% See also: transdiffreg.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2015-2016 University of Oxford
% Version: 0.4.4
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
switch (opt.Transform)
    
    case {'TranslationTransform', 'EulerTransform', 'AffineTransform'}

        narginchk(2, 2);
        
    case {'BSplineTransform'}
        
        narginchk(2, 3);
        
    otherwise
        
        error(['Transform ' opt.Transform ' not implemented'])
        
end
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
    opt.Epsilon = 0.0;
end
if (~isfield(opt, 'MaxIter'))
    opt.MaxIter = 50;
end
if (~isfield(opt, 'Alpha'))
    opt.Alpha = 0.45;
end
if (strcmp(opt.Transform, 'BSplineTransform'))
    if (~isfield(opt, 'tbsp'))
        opt.tbsp = [];
    end
end

if (opt.Alpha < 0 || opt.Alpha > 0.5)
    warning('alpha must be in the interval [0.0, 0.5]. alpha=0.5 can produce oscillations. We recommend alpha<=0.45')
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

    % stopping criterion: B-spline injectivity
    if (strcmp(opt.Transform, 'BSplineTransform'))
        
        % find control points that don't guarantee injectivity and
        % constrain them so that they do
        [info.NConstrained(I), aux] = ...
            check_bspline_injectivity_choi2000(tout + t, opt.tbsp);
        
        % apply the control point modifications to the B-spline update
        t = aux - tout;
        
    end
    
    % compose with the previous transforms to obtain a total transform
    tout = tout + t;
    
    % update neighbour transforms ...
    
    % ... of internal images
    tp(1:end-1, :) = t(2:end, :) + tp(1:end-1, :) - t(1:end-1, :);
    tm(2:end, :) = t(1:end-1, :) + tm(2:end, :) - t(2:end, :);
    
    %... of extreme images
    tp(end, :) = tm(end, :);
    tm(1, :) = tp(1, :);
    
    % stopping criterion: absolute value of each translation component
    tabs = abs(t(:));
    info.MaxAbs(I) = max(tabs);
    info.MeanAbs(I) = mean(tabs);
    info.MedAbs(I) = median(tabs);
    if (info.MaxAbs(I) <= opt.Epsilon)
        break
    end
    
    % stopping criterion: maximum number of iterations
    if (I == opt.MaxIter)
        break
    end
    
end

info.NumIter = I;

% DEBUG:
% isInjective = check_bspline_injectivity_global(tout + t, opt.tbsp);
% nnz(~isInjective)

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
    info.MedTNorm(I) = median(tnorm);
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

%% check_bspline_injectivity_global:
%
% Check B-spline global injectivity triangulating the grid and looking for
% flipped triangles.
%
% tout: (N, K)-matrix, N number of slices, K number of B-spline control
%       points
% tbsp: elastix struct with the B-spline info (we need the grid spacing)
function isInjective = check_bspline_injectivity_global(tout, tbsp)

DEBUG = false;

% number of slices
N = size(tout, 1);

% grid of control points for 0 deformation
tbsp.TransformParameters(:) = 0;
[gx0, gy0] = elastix_bspline_grid(tbsp);
tri = delaunay(gx0(:), gy0(:));

% area of the triangles
atri = trifacet_signed_area(tri, [gx0(:), gy0(:)]);

% flip triangles with negative area, so that all are positive
idx = atri < 0;
tri(idx, :) = tri(idx, end:-1:1);

% DEBUG: plot control points
if (DEBUG)
    trimesh(tri, gx0(:), gy0(:))
end

isInjective = true(1, N);
for J = 1:N

    % transfer coefficients to the B-spline
    tbsp.TransformParameters = tout(J, :);
    
    % apply transformation to the control points (note:
    aux = transformix_pts(tbsp, [gx0(:) gy0(:)]);
    gx1 = reshape(aux(:, 1), size(gx0, 1), size(gx0, 2));
    gy1 = reshape(aux(:, 2), size(gy0, 1), size(gy0, 2));
    
    % DEBUG: plot control points
    if (DEBUG)
        trimesh(tri, gx1(:), gy1(:))
    end
    
    % area of the triangles
    atri = trifacet_signed_area(tri, [gx1(:), gy1(:)]);
    
    % if any triangle is flipped, the B-spline is not injective
    isInjective(J) = all(atri > 0);
    
end

end

%% check_bspline_injectivity_choi2000:
%
% Check B-spline injectivity using sufficient conditions in Choi and Lee
% (2000). Constrain control points to guarantee injectivity.
%
% This function checks Choi and Lee (2000) conditions for injectivity of
% the B-spline. The conditions are sufficient, i.e. if any is fulfilled,
% then we know the B-spline is injective. If not, the B-spline may or may
% not be injective, but we'll assume it's not.
%
% tout: (N, K)-matrix, N number of slices, K number of B-spline control
%       points
% tbsp: elastix struct with the B-spline info (we need the grid spacing)
%
% nConstrained: number of control points that had to be constrained to
%               guarantee injectivity
function [nConstrained, tout] = check_bspline_injectivity_choi2000(tout, tbsp)

% constants from Choi and Lee (2000), to check violations
K2 = 2.046392675;
A2 = sqrt((3/2)^2 + (K2 - 3/2)^2);

% number of slices
N = size(tout, 1);

% Note: the grid has row->x, col->y coordinates, instead of the Matlab
% convention, which is the opposite

% number of B-spline control points
L = length(tbsp.TransformParameters)/2;
Nx = tbsp.GridSize(1);
Ny = tbsp.GridSize(2);

% loop slices
isInjective = false(Nx, Ny, N);
for J = 1:N
    
    %% find control points that are guaranteed to not cause injectivity 
    %% problems
    
    % transfer coefficients to the B-spline
    tbsp.TransformParameters = tout(J, :);
    
    % easy nomenclature for x, y coordinates of control points
    cx = tbsp.TransformParameters(1:L);
    cy = tbsp.TransformParameters(L+1:end);
    
    % reshape coefficients into grid
    cx = reshape(cx, Nx, Ny);
    cy = reshape(cy, Nx, Ny);
    
    % normalize control point displacement to grid spacing = 1
    cx = cx / tbsp.GridSpacing(1);
    cy = cy / tbsp.GridSpacing(2);
    
    % Theorem 1 from Choi and Lee (2000): first sufficient condition
    % for injectivity
    cond1 = (abs(cx) < 1 / K2) & (abs(cy) < (1 / K2));
    
    % Theorem 2 from Choi and Lee (2000): second sufficient condition
    % for injectivity
    cond2 = cx.^2 + cy.^2 < (1 / A2)^2;
    
    % if either condition is fulfilled, then the B-spline is injective
    isInjective(:, :, J) = cond1 | cond2;
    
    %% correct control points that potential can cause non-injectivity
    
    % find problematic points (linear index easier than row/column
    % index)
    idx = find(~isInjective(:, :, J));
    
    % correction factor for the problematic coefficients
    S = 1 ./ sqrt(cx(idx).^2 + cy(idx).^2) / (A2 - 0.01);
    
    % correct problematic control points
    cx(idx) = cx(idx) .* S;
    cy(idx) = cy(idx) .* S;
    
    % transfer corrected coefficients to the output
    aux = tout(J, :);
    aux(idx) = tbsp.GridSpacing(1)*cx(idx);
    aux(L + idx) = tbsp.GridSpacing(2)*cy(idx);
    tout(J, :) = aux;
    
end

% count number of control points that had to be constrained to guarantee
% they are injective
nConstrained = nnz(~isInjective);

end

%% check_bspline_injectivity_chun2009:
%
% Check B-spline injectivity using sufficient conditions in Chun and
% Fessler (2009).
%
% This function checks Chun and Fessler (2009) conditions for injectivity
% of the B-spline. If all the conditions are fulfilled, then we know the
% B-spline is injective. If not, the B-spline may or may not be injective,
% but we'll assume it's not.
%
% Because the conditions check grid edge lengths, it's not trivial to
% associate the conditions to vertices.
%
% tout: (N, K)-matrix, N number of slices, K number of B-spline control
%       points
% tbsp: elastix struct with the B-spline info (we need the grid spacing)
function isInjective = check_bspline_injectivity_chun2009(tout, tbsp)

DEBUG = false;

% constants from Chun and Fessler (2009), to check for violations
kx = 0.5 - 0.01;
ky = 0.5 - 0.01;

% number of slices
N = size(tout, 1);

% number of B-spline coefficients
L = length(tbsp.TransformParameters);
Nx = tbsp.GridSize(1);
Ny = tbsp.GridSize(2);

% DEBUG: plot control points
if (DEBUG)
    % grid of control points for 0 deformation
    tbsp.TransformParameters(:) = 0;
    [gx0, gy0] = elastix_bspline_grid(tbsp);
    tri = delaunay(gx0(:), gy0(:));
    
    trimesh(tri, gx0(:), gy0(:))
end

% loop slices
isInjective = true;
for J = 1:N
    
    % transfer coefficients to the B-spline
    tbsp.TransformParameters = tout(J, :);
    
    % easy nomenclature for x, y coordinates of control points
    cx = tbsp.TransformParameters(1:L/2);
    cy = tbsp.TransformParameters(L/2+1:end);
    
    % reshape coefficients into grid
    cx = reshape(cx, Nx, Ny);
    cy = reshape(cy, Nx, Ny);
    
    % transpose grid so that we have x->cols, y->rows
    cx = cx';
    cy = cy';
    
    % grid spacing
    mx = tbsp.GridSpacing(1);
    my = tbsp.GridSpacing(2);
    
    % conditions (i->x, j->y):
    
    % -mx * kx <= c_{i+1,j}^x - c_{i,j}^x
    % we also add a column to the right to account for the lost column
    cond = (diff(cx, 1, 2) >= -mx * kx);
    isInjective = isInjective && all(cond(:));
    
    % -my * ky <= c_{i,j+1}^y - c_{i,j}^y
    % we also add a row to the bottom to account for the lost row
    cond = (diff(cy, 1, 1) >= -my * ky);
    isInjective = isInjective && all(cond(:));
    
    % |c_{i,j+1}^x - c_{i,j}^x| <= mx * kx
    % we also add a row to the bottom to account for the lost row
    cond = (abs(diff(cx, 1, 1)) <= mx * kx);
    isInjective = isInjective && all(cond(:));
    
    % |c_{i+1,j}^y - c_{i,j}^y| <= my * ky
    % we also add a column to the right to account for the lost column
    cond = (abs(diff(cy, 1, 2)) <= my * ky);
    isInjective = isInjective && all(cond(:));

end

end
