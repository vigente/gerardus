function [ CF, I_grad] = CF_dictionary_3D(Image, k_space_undersampled, mask, samp_pdf, beta_TV, V, I_sparse, beta_sparse, TE_mask, T2, AT)
% returns only the DC and TV residuals

if nargout < 2
    get_gradient = false;
else
    get_gradient = true;
end 

sz = size(Image);

% Sparsity term - In CS reconstruction, we could let the sparsity residual
% (Image - I_sparse) be represented in the absolutely least sparse way -
% imagine the dictionary had extra columns with a 1 for each diffusion
% direction, i.e. [1 0 0 0 ...; 0 1 0 0 ...; 0 0 1 0 ...]' 
Sparse_grad = abs(Image) - abs(I_sparse);

Sparse_value = full(sum(abs(V(:)))) + ... % L1 norm of sparse signal
                sum(abs(Sparse_grad(:))); % L1 norm of errors
        

% Total variation
% TV applies only to the magnitudue image
[TV_value, TV_grad] = forward_TV_DWI(abs(Image));

  
% Data consistency

% % get k-space
% %disp('Getting k-space')
% k_space = fft3c(Image); % performs fft only in first 3 dims

% get k-space accounting for T2-weighting
k_space = zeros(size(Image));
TE = unique(TE_mask);
TE = setdiff(TE, 0);
for i = 1:length(TE)
    
    % which parts of k-space are needed? (TE_mask is 4D)
    TE_mask_k = TE_mask == TE(i);

    % get an idea of what the decayed image should be
    I_decayed = bsxfun(@times, Image, exp(- (TE(i) - TE(1)) / T2));

    % get the fourier transform of the decayed image (a portion of it
    % will be used later)
    K_decayed = fft3c(I_decayed);
        
    % put the right values into k-space
    k_space(TE_mask_k) = K_decayed(TE_mask_k);
    
end



% get the data consistency
%[DC_value, DC_grad] = data_consistency(k_space, k_space_undersampled, mask, samp_pdf);
[DC_value, DC_grad] = data_consistency(k_space, k_space_undersampled, mask); % changed for T2 stuff



% Add the errors up
CF = {};
CF.CF = DC_value + beta_TV*TV_value + beta_sparse * Sparse_value;
CF.DC = DC_value;
CF.TV = beta_TV*TV_value;
CF.L1 = beta_sparse * Sparse_value;


% Get the gradient if needed


if get_gradient
    
    %disp('Getting TV residuals')
    % TV residuals
    TV_resids = inverse_TV_DWI(TV_grad ./ (abs(TV_grad) + 1E-15));
    %TV_resids = inverse_TV_DWI(TV_grad / 6 + 1E-15);
    
    
    % Give the TV resids their proper phase
    I_phase = angle(Image);
    TV_resids = TV_resids .* exp(1i * I_phase);
    
    clear TV_grad
    
    
    % DC residuals
    %DC_grad_I = ifft3c(DC_grad) * 2;
    
    DC_grad_I = zeros(size(Image));   
    for i = 1:length(TE)
        
        % which parts of k-space are needed? (TE_mask is 4D)
        TE_mask_k = TE_mask == TE(i);

        % convert back to image
        I_resids = ifft3c(DC_grad .* TE_mask_k);

        I_resids = bsxfun(@rdivide, I_resids, exp(- (TE(i)-TE(1)) / T2));
        
        DC_grad_I = DC_grad_I + I_resids;
    end
    DC_grad_I = DC_grad_I * 2;
    
    
           
            
    Sparse_grad = inverse_basis_function(spones(V), AT, sz) + ...
            Sparse_grad ./ (abs(Sparse_grad) + 1E-15);


    
    
    % give the sparse residuals their correct phase
    Sparse_grad = Sparse_grad .* exp(1i * I_phase);
    
    
    % Add it up
    I_grad = DC_grad_I + beta_TV * TV_resids + beta_sparse * Sparse_grad;
    
    
    figure(200);
    subplot(4,1,1); hold off;
    plot(real(squeeze(DC_grad_I(100, 90, 60, :))), 'b'); 
    hold on; plot(imag(squeeze(DC_grad_I(100, 90, 60, :))), 'g'); 
    hold on; plot(abs(squeeze(DC_grad_I(100, 90, 60, :))), 'r'); 
    title('DC')
    subplot(4,1,2); hold off;
    plot(real(squeeze(beta_TV*TV_resids(100, 90, 60, :))), 'b'); 
    hold on; plot(imag(squeeze(beta_TV*TV_resids(100,90, 60, :))), 'g'); 
    hold on; plot(abs(squeeze(beta_TV*TV_resids(100,90, 60, :))), 'r'); 
    title('TV')
    subplot(4,1,3); hold off;
    plot(real(squeeze(beta_sparse * Sparse_grad(100,90, 60, :))), 'b'); 
    hold on; plot(imag(squeeze(beta_sparse * Sparse_grad(100, 90, 60, :))), 'g'); 
    hold on; plot(abs(squeeze(beta_sparse * Sparse_grad(100, 90, 60, :))), 'r'); 
    title('L1')
    subplot(4,1,4); hold off;
    plot(real(squeeze(I_grad(100,60, 60, :))), 'b'); 
    hold on; plot(imag(squeeze(I_grad(100, 90, 60, :))), 'g'); 
    hold on; plot(abs(squeeze(I_grad(100, 90, 60, :))), 'r'); 
    title('Sum')
    pause(1)
    
end
