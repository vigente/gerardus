function [ V ] = forward_basis_function(Image_space, AT, sz, noise_level)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

A_L1 = AT';
 l_norm = 2;
 A_L2 = A_L1 ./ repmat((sum(A_L1.^l_norm ,1).^(1/l_norm)), [size(A_L1,1), 1]);

AT = A_L2';
A = double(AT');
 
S = reshape(Image_space, [prod(sz(1:end-1)), sz(end)]);

S_hat = double(S');

param.mode=1; % 2
param.pos=1; % 1
param.lambda= noise_level * sz(end); %0.005; 
param.cholesky = 0;
% tic
% V_orig = mexLasso(S_hat,double(AT)',param)';
% toc

%---- re-weighting in parfor

s = matlabpool('size');
if s == 0
    matlabpool open
end

batchsize = 2560;
dictsize = size(AT,1);
sigsize = size(S,2);
numbatches = floor(size(S,1) / batchsize);

bsnb = batchsize * numbatches; % = num voxels
bsds = batchsize * dictsize; % = size after reshaping

V_batch = sparse(numbatches, bsds);
V_batch_L1 = sparse(numbatches, bsds);

S_batch = reshape(S_hat(:, 1:bsnb), [batchsize*sigsize, numbatches])';

parfor i = 1:numbatches
    
    S_batch_reshaped = reshape(S_batch(i,:), [sigsize, batchsize]);
    
    Y = mexLasso(S_batch_reshaped, A, param);
    
    V_batch_L1(i,:) = reshape(Y, [1, bsds]);
    
    % weighting once...
    W = 1./(Y + 0.01); % this may need lots of ram
    Y = mexLassoWeighted(S_batch_reshaped, A, W, param);
        
    % weighting twice
    W = 1./(Y + 0.01);
    Y = mexLassoWeighted(S_batch_reshaped, A, W, param);
    
    V_batch(i,:) = reshape(Y, [1, bsds]);

end

V_L1 = reshape(V_batch_L1', [dictsize, bsnb])';
V = reshape(V_batch', [dictsize, bsnb])';

% find the voxels that failed the reweighting step
[r, ~] = find(V < 0);
r = unique(r);
if ~isempty(r)
    disp([num2str(length(r)/length(V)*100) '% of voxels failed the reweighting step'])
end

V(r, :) = V_L1(r, :);

I_L1 = inverse_basis_function(V_L1, AT, [prod(sz(1:3)), sz(4)]);
I_L0 = inverse_basis_function(V, AT, [prod(sz(1:3)), sz(4)]);


cf_L1 = full(sum(abs(V_L1),2)) + ... % L1 norm of sparse signal
                sum(I_L1 - S,2); % L1 norm of errors

cf_L0 = full(sum(abs(V),2)) + ... % L1 norm of sparse signal
    sum(I_L0 - S,2); % L1 norm of errors

r = cf_L1 < cf_L0;

if ~isempty(r)
    disp(['A further ' num2str(sum(r)/length(V)*100) '% of voxels were better before reweighting'])
end

V(r, :) = V_L1(r, :);

V(isnan(V)) = 0;

