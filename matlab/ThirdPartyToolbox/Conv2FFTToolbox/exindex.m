function arr = exindex(arr, varargin)
%EXINDEX extended array indexing
%   ARROUT = EXINDEX(ARRIN, S1, S2, ...) indexes a virtual array made by
%   extending ARRIN with zeros in all directions, using subscripts S1, S2
%   etc.
%
%   ARROUT = EXINDEX(ARRIN, S1, R1, S2, R2, ...) extends ARRIN using rule
%   R1 on the first dimension, R2 on the second dimension etc.
%
%   ARROUT = EXINDEX(ARRIN, S1, S2, ..., R) extends ARRIN using rule R on
%   every dimension.
%
%   Subscripts
%   ----------
%
%   Broadly, if V is the virtual extended array, ARROUT = V(S1, S2, ...)
%
%   The elements of the subscript arguments S1, S2 etc must be integers.
%   They need not be positive and are not restricted in any way by the size
%   of ARRIN. Logical indexing and linear indexing are not supported.
%
%   There must be at least one subscript argument for each dimension of
%   ARRIN as reported by NDIMS, except that row and column vectors may have
%   1 or 2 subscripts. A single subscript is taken to refer to the
%   dimension along which the vector lies, as in normal vector indexing.
%   Scalars require 2 subscripts. If there are more subscripts than
%   dimensions, ARRIN is taken to have trailing singleton dimensions, as in
%   normal array indexing.
%
%   The number of dimensions of ARROUT will be the number of subscript
%   arguments, though trailing singleton dimensions will, as usual, be
%   suppressed. The size of ARROUT is given by the normal Matlab rules for
%   the result of indexing into ARRIN: that is
%
%       size(ARROUT) = size( ARRIN(ones(size(S1)), ones(size(S2)), ...) )
%
%   A subscript argument may be the string ':'. This behaves like a colon
%   in ordinary subscripting: a colon for the K'th subscript stands for
%   1:size(ARRIN, K). The 'end' keyword is not supported.
%
%   Rules
%   -----
%
%   Each rule may be one of the following:
%
%   A scalar cell: ARRIN is padded with elements equal to the contents of
%   the cell. The class of the cell contents must be compatible with the
%   class of ARRIN.
%
%       If different constants are used on different dimensions, padding is
%       done in the order of the subscripts. For example, a 2D array is
%       extended first in the row index direction and then in the column
%       index direction. For all other cases, the order in which dimensions
%       are extended has no effect.
%
%   'circular': ARRIN is extended with copies of itself; i.e. V is tiled
%   with ARRIN.
%
%   'symmetric': ARRIN is extended with copies of itself with reflection at
%   its boundaries; i.e. V is tiled with [ARRIN fliplr(ARRIN);
%   flipud(ARRIN) fliplr(flipud(ARRIN))].
%
%   'replicate': ARRIN is extended by copying its border elements; i.e. an
%   element of V is equal to the nearest element of ARRIN.
%
%   If no rule is given, padding is with zeros.
%
%   Examples
%   --------
%
%   Pad a 2D matrix with K extra rows and columns with reflection on both
%   axes:
%
%       b = exindex(a, 1-k:size(a,1)+k, 1-k:size(a,2)+k, 'symmetric');
%
%   Circularly shift a 2D matrix by R rows downwards and C columns
%   rightwards:
%
%       b = exindex(a, 1-r:size(a,1)-r, 1-c:size(a,2)-c, 'circular');
%
%   Force a row or column vector to be 1024 elements long, trimming or
%   padding with zeros as necessary:
%
%       u = exindex(v, 1:1024);
%
%   The same, with a non-zero padding value:
%
%       u = exindex(v, 1:1024, {-1});   % note constant in cell
%
%   Truncate or extend all the rows of a matrix to 1024 columns:
%
%       b = exindex(a, ':', 1:1024);
%
%   Extend a 2-D array into the third dimension by copying it:
%
%       b = exindex(a, ':', ':', 1:3, 'replicate');
%
%   Pad a 1-D cell array with cells containing the empty matrix:
%
%       cellout = exindex(cellin, 0:10, {{[]}}); 
%
%   See also: padarray, circshift, repmat

% Copyright David Young 2010

% Sort out arguments
[exindices, rules, nd, sz] = getinputs(arr, varargin{:});
consts = cellfun(@iscell, rules);  % Check for constants, as can be
constused = any(consts);           % more efficient if there are none

% Setup for constant padding
if constused
    tofill = cell(1, nd);
end

% Main loop over subscript arguments, transforming them into valid
% subscripts into arr using the rule for each dimension
if constused
    for i = 1:nd
        [exindices{i}, tofill{i}] = extend(exindices{i}, rules{i}, sz(i));
    end
else % no need for information for doing constants
    for i = 1:nd
        exindices{i} = extend(exindices{i}, rules{i}, sz(i));
    end
end

% Create the new array by indexing into arr. If there are no constants,
% this does the whole job
arr = arr(exindices{:});

% Fill areas that need constants
if constused
    % Get full range of output array indices
    ranges = arrayfun(@(x) {1:x}, size(arr));
    for i = nd:-1:1    % order matters
        if consts(i)
            ranges{i} = tofill{i};      % don't overwrite original
            c = rules{i};               % get constant and fill ...
            arr(ranges{:}) = c{1};      % we've checked c is scalar
            ranges{i} = ~tofill{i};     % don't overwrite
        end
    end
end

end

% -------------------------------------------------------------------------

function [exindices, rules, nd, sz] = getinputs(arr, varargin)
% Sort out and check arguments. Inputs are as given in the help comments
% for exindex. Outputs are cell arrays; each element of exindices is a
% set of integer extended indices which has been checked for validity; each
% element of rules is a rule which has not been checked for validity.

% Use index/rules arguments only to establish no. dimensions - ndims(arr)
% is no use, as trailing singleton dimensions truncated and vectors can be
% 2D or 1D
nd = length(varargin);
if nd == 0
    error('exindex:missingargs', 'Not enough arguments');
elseif nd == 1
    exindices = varargin;
    rules = {{0}};
elseif ~(isnumeric(varargin{2}) || strcmp(varargin{2}, ':'))
    % have alternating indices and rule
    nd = nd/2;
    if round(nd) ~= nd
        error('exindex:badnumargs', ...
            'Odd number of arguments after initial index/rule pair');
    end
    exindices = varargin(1:2:end);
    rules = varargin(2:2:end);
elseif nd > 2 && ~(isnumeric(varargin{end}) || strcmp(varargin{end}, ':'))
    % have a general rule at end
    nd = nd - 1;
    exindices = varargin(1:nd);
    [rules{1:nd}] = deal(varargin{end});
else
    % no rule is specified
    exindices = varargin;
    [rules{1:nd}] = deal({0});
end

% Sort out mismatch of apparent array size and number of dimensions
% indexed
sz = size(arr);
ndarr = ndims(arr);
if nd < ndarr
    if nd == 1 && ndarr == 2
        % Matlab allows vectors to be indexed with a single subscript and
        % to retain their shape. In all other cases (including scalars) a
        % single subscript causes the output to take the same shape as the
        % subscript array - we can't deal with this.
        if sz(1) == 1 && sz(2) > 1
            % have a row vector
            exindices = [{1} exindices {1}];
            rules = [rules rules];  % 1st rule doesn't matter
        elseif sz(2) == 1 && sz(1) > 1
            % have a column vector
            exindices = [exindices {1}];
            rules = [rules rules];  % 2nd rule doesn't matter
        else
            error('exindex:wantvector', ...
                'Only one index but array is not a vector');
        end
    else
        error('exindex:toofewindices', ...
            'Array has more dimensions than there are index arguments');
    end
    nd = 2;
elseif nd > ndarr
    % Effective array size
    sz = [sz ones(1, nd-ndarr)];
end

% Expand any colons now to simplify checking.
% It's tempting to allow the 'end' keyword here: easy to substitute the
% size of the dimension. However, to be worthwhile it would be necessary to
% use evalin('caller',...) so that expressions using end could be given as
% in normal indexing. This would mean moving the code up to exindex itself,
% and evalin makes for inefficiency and fragility, so this hasn't been
% done.
colons = strcmp(exindices, ':');
if any(colons)  % saves a little time
    exindices(colons) = arrayfun(@(x) {1:x}, sz(colons));
end

% Check the indices (rules are checked as required in extend)
checkindex = @(ind) validateattributes(ind, {'numeric'}, ...
    {'integer'}, 'exindex', 'index');
cellfun(checkindex, exindices);

end

% -------------------------------------------------------------------------

function [ind, tofill] = extend(ind, rule, s)
% The core function: maps extended array subscripts into valid input array
% subscripts.

if ischar(rule)    % pad with rule
    
    tofill = [];  % never used
    switch rule
        case 'replicate'
            ind = min( max(1,ind), s );
        case 'circular'
            ind = mod(ind-1, s) + 1;
        case 'symmetric'
            ind = mod(ind-1, 2*s) + 1;
            ott = ind > s;
            ind(ott) = 2*s + 1 - ind(ott);
        otherwise
            error('exindex:badopt', 'Unknown option');
    end
    
elseif iscell(rule) && isscalar(rule)     % pad with constant
    
    % The main messiness is due to constant padding. This can't be done
    % with indexing into the original array, but we want the indexing
    % structure to be preserved, so for now we index to element 1 on each
    % dimension, and record the indices of the regions that need to be
    % fixed.
    
    tofill = ind < 1 | ind > s;
    ind(tofill) = 1;
    
else
    
    error('exindex:badconst', 'Expecting string or scalar cell');
    
end

end


