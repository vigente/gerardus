function res = ifftNc(x)
% circular ifft

res = sqrt(length(x(:)))*ifftshift(ifftn(fftshift(x)));

