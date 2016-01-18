function K = fft3c(x)
% K = FFT3C(X) performs a circular fft in the first 3 dimensions

% X can have any dimensionality.
% Matlab's fft assumes DC is at index 1, so this function should be used 
% when you have DC at end/2 + 1



% Note: this function seems to not work so well with odd sized images

% if the size is 3, use fftNc - it is much faster 
if length(size(x)) == 3
    K = fftNc(x);
    return
end


sz = size(x);
if length(sz) < 3
    sz(3) = 1;
end

% this approach shifts all dimensions (even the ones you don't want
% shifted), performs the ffts, and then shifts back. It is faster than
% shifting the three dimensions that you do want individually.

K = 1/sqrt(sz(1)*sz(2)*sz(3)) * fftshift(fft(fft(fft(ifftshift(x), [], 3), [], 2), [], 1));


% %res = 1/sqrt(length(x(:)))*fftshift(fft2(ifftshift(x)));
% res = 1/sqrt(length(x(:)))*fftshift_complex( ...
%         fft2(fft(...
%         ifftshift_complex(x)...
%         , [], 3))...
%     );
% 
