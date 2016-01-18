function [ DC, grad_DC] = data_consistency(k_recon, k, mask, pdf)
% returns the data consistency (similarity to sampled k_space) of k_recon

% where mask = 0, the gradient is also zero since this should not influence
% the gradient. otherwise, the gradient is simply the gradient of DC with respect to k(i)

if nargin < 4
    pdf = ones(size(mask));
end

k_diff = (k_recon - k) .* mask ./ pdf;

DC = k_diff(:)' * k_diff(:);

grad_DC = single(k_diff);



% % x is the current value
% x = k_recon(mask);
% 
% % a is the ground truth
% a = k(mask);
% 
% % b is everything except the i-th term
% b = sum_of_squares - abs(x-a).^2;
% 
% % derivative of sqrt((x - a)^2 + b)
% grad_DC = zeros(size(k));
% grad_DC(mask) = - (a - x) ./ (sqrt((abs(a - x)).^2 + b));
% 
% % anywhere that is nan is zero
% grad_DC(isnan(grad_DC)) = 0;
