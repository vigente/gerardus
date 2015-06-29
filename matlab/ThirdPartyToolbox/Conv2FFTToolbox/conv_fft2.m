function y = conv_fft2(x, m, shape)
%CONV_FFT2 Two dimensional convolution via the FFT.
%   Y = CONV_FFT2(X, M) performs the 2-D convolution of matrices X and M.
%   If [mx,nx] = size(X) and [mm,nm] = size(M), then size(Y) =
%   [mx+mm-1,nx+nm-1]. Values near the boundaries of the output array are
%   calculated as if X was surrounded by a border of zero values.
%
%   Y = CONV_FFT2(X, M, SHAPE) where SHAPE is a string returns a
%   subsection of the 2-D convolution with size specified by SHAPE:
%
%       'full'    - (default) returns the full 2-D convolution,
%       'same'    - returns the central part of the convolution
%                   that is the same size as X (using zero padding),
%       'valid'   - returns only those parts of the convolution
%                   that are computed without the zero-padded
%                   edges, size(Y) = [mx-mm+1,nx-nm+1] when
%                   size(X) > size(M),
%       'wrap'    - as for 'same' except that instead of using
%                   zero-padding the input X is taken to wrap round as
%                   on a toroid.
%       'reflect' - as for 'same' except that instead of using
%                   zero-padding the input X is taken to be reflected
%                   at its boundaries.
%
%   For shape options 'full', 'same' and 'valid' this function should
%   return the same result as CONV2, to within a small tolerance. For all
%   shape options, this function should return the same result as CONVOLVE2
%   (available via the MATLAB Central File Exchange), with no TOL argument,
%   to within a small tolerance.
%
%   CONV_FFT2 uses multiplication in the frequency domain to compute the
%   convolution. It may be faster than CONV2 and CONVOLVE2 for masks above
%   a certain size. This should be checked experimentally for any given
%   application and system.
%
%   See also CONV2, CONVOLVE2, FILTER2.

% Copyright David Young, April 2011

error(nargchk(2,3,nargin,'struct'));
if nargin < 3
    shape = 'full';    % shape default as for CONV2
end;

[x, m, fsize] = padarrays(x, m, shape);

% no need to trap case of real x and m - fft2 handles efficiently
y = ifft2(fft2(x) .* fft2(m));   % central operation, basic form

% trim to correct output size
if ~isequal(fsize, size(y))
    y = y(1:fsize(1), 1:fsize(2));
end

end

function [x, m, fsize] = padarrays(x, m, shape)
% Pad arrays to make them the same size and allow for boundary effects

xsize = size(x);
msize = size(m);

switch shape
    
    case 'wrap'
        fsize = xsize;
        % ensure x no smaller than m
        if any(msize > xsize)  && ~isempty(x)
            x = repmat(x, ceil(msize ./ size(x)));
            xsize = size(x);
        end
        % pad m with zeros
        if any(msize < xsize)  % test, as user may have optimised already
            m = exindex(m, 1:xsize(1), 1:xsize(2), {0});
        end
        % recentre m so that y(1,1) corresponds to mask centred on x(1,1)
        mc = 1 + floor(msize/2);
        me = mc + xsize - 1;
        m = exindex(m, mc(1):me(1), mc(2):me(2), 'circular');
    
    case 'full'
        fsize = xsize + msize - 1;  % enough room for no overlap
        x = exindex(x, 1:fsize(1), 1:fsize(2), {0});
        m = exindex(m, 1:fsize(1), 1:fsize(2), {0});
    
    case 'valid'
        fsize = xsize - msize + 1;
        % pad m with zeros (don't test first, as likely to be needed)
        m = exindex(m, 1:xsize(1), 1:xsize(2), {0});
        % shift m so that y(1,1) corresponds to mask just inside x
        me = msize + xsize - 1;
        m = exindex(m, msize(1):me(1), msize(2):me(2), 'circular');

    case 'same'
        fsize = xsize;
        mmid = floor(msize/2);
        xsize = xsize + mmid;   % border to avoid edge effects
        x = exindex(x, 1:xsize(1), 1:xsize(2), {0});
        m = exindex(m, 1:xsize(1), 1:xsize(2), {0});
        % recentre m so that y(1,1) corresponds to mask centred on x(1,1)
        mc = 1 + mmid;
        me = mc + xsize - 1;
        m = exindex(m, mc(1):me(1), mc(2):me(2), 'circular');
        
    case 'reflect'
        fsize = xsize;
        xsize = xsize + msize - 1;   % border to avoid edge effects
        xc = 1 - floor((msize-1)/2);
        xe = xc + xsize - 1;
        x = exindex(x, xc(1):xe(1), xc(2):xe(2), 'symmetric');
        m = exindex(m, 1:xsize(1), 1:xsize(2), {0});
        % recentre m so that y(1,1) corresponds to mask centred on x(1,1)
        me = msize + xsize - 1;
        m = exindex(m, msize(1):me(1), msize(2):me(2), 'circular');

    otherwise
        error('conv_fft2:badshapeopt', 'Unrecognised shape option: %s', shape);
end
end

