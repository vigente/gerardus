function res = fftNc(x)
% circular fft

res = 1/sqrt(length(x(:)))*fftshift(fftn(ifftshift(x)));
