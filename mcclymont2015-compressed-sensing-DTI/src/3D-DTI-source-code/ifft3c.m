function res = ifft3c(x)
% performs the fft in the first 3 dimensions only
% if your data set is 3D anyway, it is much faster to use ifftNc


if length(size(x)) == 3
    res = ifftNc(x);
    return
end
% 
% 
% % to do: find a faster way of doing this
% res = zeros(size(x));
% for i = 1:size(x,4)
%     res(:,:,:,i) = ifftNc(x(:,:,:,i));
% end

sz = size(x);
res = sqrt(sz(1)*sz(2)*sz(3)) * ifftshift(ifft(ifft(ifft(fftshift(x), [], 3), [], 2), [], 1));