% Generate 2D phantom
close all
clear all
clc

[X, Y] = ndgrid(1:160, 1:160);

center_x = 80;
center_y = 80;

Buffer = ((X-center_x).^2 + (Y-center_y).^2) < 20^2;

Tissue = ((X-center_x).^2 + (Y-center_y).^2) < 50^2;
Tissue = Tissue & ~Buffer;

Gel = ((X-center_x).^2 + (Y-center_y).^2) < 70^2;
Gel = Gel & ~Tissue & ~Buffer;


% define T2 weighting
T2 = zeros(size(X)) + inf;
T2(Buffer) = 0.04;
T2(Tissue) = 0.024;
T2(Gel) = 0.03;

TE = [ 0.0093    0.0142    0.0191    0.0240    0.0289    0.0338    0.0387    0.0437];

% get the b values

bval = procpar_to_bval(['REDACTED']); % this is a 3*3*N array
% re-order the sequence
dwi_order = [1, 9, 17, 25, setdiff(1:34, [1, 9, 17, 25])]; % just re-ordering to put b=0 first
bval = bval(:,:,dwi_order);

Bv=squeeze([bval(1,1,:),2*bval(1,2,:),2*bval(1,3,:),bval(2,2,:),2*bval(2,3,:),bval(3,3,:)])';


I_DTI = zeros([size(X), 34]);


for i = 1:size(X,1)
    for j = 1:size(X,2)

        if Tissue(i,j)
        
            % radial vector (2D)
            v_r = [i - center_x, j - center_x, 0];

            v_circ = cross(v_r, [0 0 1]);
            
            % normalise
            v_circ = v_circ / norm(v_circ);

            % v_circ needs to be given a z component proportional to its
            % distance
            d = sqrt(sum(v_r.^2));
            
            ha = d * 6 - 208; % desired helix angle
            %ha = d * 6 - 210; % desired helix angle
            
            dp = cos(pi/2 - ha/180*pi); % required dot product between v_cl and v_l
            
            % kx = v_circ(1), ky = v_circ(2), (kx)^2 + (ky)^2 = 1 - dp^2
            ksq = (1 - dp^2) / (v_circ(1)^2 + v_circ(2)^2);
            k = sqrt(ksq);
            
            
            v_primary = [v_circ(1)*k, v_circ(2)*k, dp];
            
            v_tertiary = v_r / norm(v_r);

            v_secondary = cross(v_primary, v_tertiary);
            
            v_primary = v_primary([3 1 2]);
            v_secondary = v_secondary([3 1 2]);
            v_tertiary = v_tertiary([3 1 2]);


            ev1 = 1.3E-3;
            ev2 = 1E-3;
            ev3 = 0.7E-3;

            ev = diag([ev1, ev2, ev3]);
            evect = [v_primary(:), v_secondary(:), v_tertiary(:)];

            % get diffusion tensor
            D = evect * ev * pinv(evect);


            % convert to a vector (should be diagonally symmetric)
            M = [D(1,1), D(1,2), D(1,3), D(2,2), D(2,3), D(3,3)]';


            % get the signal for this tensor

            for n = 1:size(Bv,1)

                Bv_n = Bv(n,:);

                I_DTI(i,j,n) = 0.8 * exp(-Bv_n * M);

            end
        
        elseif Gel(i,j) || Buffer(i,j)
            
            v_primary = rand([1 3]);
            v_secondary = rand([1 3]);
            v_tertiary = rand([1 3]);
            v_primary = v_primary / norm(v_primary);
            v_secondary = v_secondary / norm(v_secondary);
            v_tertiary = v_tertiary / norm(v_tertiary);
            
            if Gel(i,j)
                ev1 = 2.2E-3;
                ev2 = 2.2E-3;
                ev3 = 2.2E-3;
            else
                ev1 = 2.3E-3;
                ev2 = 2.3E-3;
                ev3 = 2.3E-3;
            end
            
            % increase slightly in buffer
            
            ev = diag([ev1, ev2, ev3]);
            evect = [v_primary(:), v_secondary(:), v_tertiary(:)];

            % get diffusion tensor
            D = evect * ev * pinv(evect);


            % convert to a vector (should be diagonally symmetric)
            M = [D(1,1), D(1,2), D(1,3), D(2,2), D(2,3), D(3,3)]';


            % get the signal for this tensor

            for n = 1:size(Bv,1)

                Bv_n = Bv(n,:);

                I_DTI(i,j,n) = exp(-Bv_n * M);

            end
            
        end
            
    end
end


SNR = 60;
I_DTI = I_DTI * SNR;

[ DT_orig, FA_orig, ADC_orig, VectorField_orig, EigVals_orig] = fit_DT( I_DTI + eps, bval);
DT_plot = zeros([3 3 size(DT_orig,1) size(DT_orig,2)]);

rr = [3 1 2];

DT_plot(rr(1),rr(1),:,:) = DT_orig(:,:,1);
DT_plot(rr(1),rr(2),:,:) = DT_orig(:,:,2);
DT_plot(rr(1),rr(3),:,:) = DT_orig(:,:,3);
DT_plot(rr(2),rr(1),:,:) = DT_orig(:,:,2);
DT_plot(rr(2),rr(2),:,:) = DT_orig(:,:,4);
DT_plot(rr(2),rr(3),:,:) = DT_orig(:,:,5);
DT_plot(rr(3),rr(1),:,:) = DT_orig(:,:,3);
DT_plot(rr(3),rr(2),:,:) = DT_orig(:,:,5);
DT_plot(rr(3),rr(3),:,:) = DT_orig(:,:,6);


rad_vec = zeros([160 160 3]);
rad_vec(:,:,1) = (X - center_x);
rad_vec(:,:,2) = (Y - center_y);
rad_vec = bsxfun(@rdivide, rad_vec, sqrt(sum(rad_vec.^2, 3)));

long_vec = zeros([160 160 3]);
long_vec(:,:,3) = 1;

circ_vec = cross(rad_vec, long_vec, 3);

HA_orig = compute_helix_angle_on_image(rad_vec, circ_vec, long_vec, VectorField_orig(:,:,[2 3 1],1));


return

% generate dictionary
param.iter = 200;
param.K = 100;
param.lambda = 1E-2;%2%34;
param.mode = 2;%4;
param.posAlpha=1;
param.posD=1;
param.clean = 1;
param.batchsize = 512;
param.verbose = 0;

sz = size(I_DTI);
S = reshape(abs(I_DTI), [prod(sz(1:end-1)), sz(end)]);


[U, model] = mexTrainDL(S(randperm(size(S,1)), :)', param);
[val, idx] = sort(U(1,:));
U = U(:,idx);
f = find(U(1,:) > 0, 1, 'first');
A = U(:, f:end);
AT = A';


for reps = 1:10
    
    % add noise to k-space

    K_noise = 1/8 * (randn(size(I_DTI)) + 1i*randn(size(I_DTI)));

    for accel_factor = 2:6

        [mask_2D, PE1, PE2, samp_pdf_dir] = generate_FSE_DTI_mask(accel_factor, 0.15, 160, 160, 8);

        TE_line = max(mask_2D,[], 2);

        % change this for prospectve vs retrospective
        a = TE_line';
        %a = imresize([8:-1:1, 1:8], [1 160], 'nearest'); 
        
        TE_mask = repmat(abs(a'), [1 160]); % retrospective one
        % don't forget to change the signal generation too

        TE_mask = repmat(TE_mask, [1 1 34]);

        for i = 1:8
            TE_mask(TE_mask == i) = TE(i);
        end



        % generate k space with appropriate T2 weighting

        K = zeros(size(I_DTI));

        for i = 1:length(TE)

    %         inds = floor((0:159) / 10);
    %         tops = inds == (8 - i);
    %         bottoms = inds == (7 + i);


            % which parts of k-space are needed? (TE_mask is 4D)
            TE_mask_k = TE_mask == TE(i);

            % get an idea of what the decayed image should be
            %I_decayed = bsxfun(@times, Image, exp(- (TE(i)) ./ T2));
            I_decayed = bsxfun(@times, I_DTI, exp(- (TE(i) - TE(1)) ./ T2));

            % get the fourier transform of the decayed image (a portion of it
            % will be used later)
            K_decayed = fft2c(I_decayed);

            % put the right values into k-space
            K(TE_mask_k) = K_decayed(TE_mask_k);

    %         tops = TE_line == i;
    %         bottoms = TE_line == i;
    % 
    %         % get an idea of what the decayed image should be
    %         I_decayed = bsxfun(@times, I_DTI, exp(- (TE(i)-TE(1)) ./ T2));%
    % 
    % 
    %         % get the fourier transform of the decayed image (a portion of it
    %         % will be used later)
    %         K_decayed = fft2c(I_decayed);
    % 
    %         % put the right values into k-space
    %         K(tops,:,:) = K_decayed(tops,:,:);
    %         K(bottoms,:,:) = K_decayed(bottoms,:,:);


        end


        K_clean = K;

        % Add noise
        K = K + K_noise; % adding noise
        % SNR ~= 60




        % undersample K space
        central_region = 0.15;

        samp_pdf = zeros(size(K));
        mask = zeros(size(K));
        for i = 1:size(K,3)
            [mask_2D, PE1, PE2, samp_pdf_dir] = generate_FSE_DTI_mask(accel_factor, central_region, 160, 160, 8);
            mask(:,:,i) = mask_2D ~= 0;
            samp_pdf(:,:,i) = samp_pdf_dir;
        end
        mask = mask == 1;
        samp_pdf(:) = 1;

        K_us = K .* mask;

        % generate TE mask
        %TE_mask = repmat(max(mask_2D, [], 2), [1 160]); % prospective one

    %     a = (floor((-80:79) / 10));
    %     a(a >= 0) = a(a >= 0) + 1;







        disp(' Not running sliding window')
        k_zero = DTI_sliding_window(K_us(:,:,1:4), bval(:,:,1:4));

        k_high = DTI_sliding_window(K_us(:,:,5:end), bval(:,:,5:end));

        I = ifft2c(cat(3, k_zero, k_high));






        clear k_high
        clear k_zero

        % for acc = 6, TV = 1E-1 works well


        beta_TV = 1;%/(accel_factor^2);%1E-1;

        beta_L1 = 10;%10/(accel_factor^2);  

        x_n_minus_one = zeros(size(I));
        s_n_minus_one = zeros(size(I));

        cg_alpha = 0.5;
        CF_old.CF = inf;

        sz = size(I);

        set_beta_zero = false;

        sparsity_error = 1;

        FA_rmse = [];
        ADC_rmse = [];
        HA_rmse = [];
        
        FA_mean = [];
        ADC_mean = [];
        FA_std = [];
        ADC_std = [];

        for itn = 1:200

            disp(['Iteration ' num2str(itn)])

            % 1. get best estimate of V
            % V* = argmin_V { ||I - D*V||_2^2 + lambda_1||V||_1 }

            if itn > 1
               sparsity_error = mean(mean(sum((abs(I) - abs(I_sparse)).^2, 3))) / sz(end) * 2;
            end

            if beta_L1 > 0

                if itn > 1
                    V_old = V;
                    I_sparse_old = I_sparse;
                end

                disp(['Running forward transform, allowable error = ' num2str(sparsity_error)])
                V = forward_basis_function(abs(I), AT, sz, sparsity_error);

                figure(5); imagesc(reshape(sum(V~=0,2), [160 160]), [0 15]); colorbar; drawnow

                I_sparse = inverse_basis_function(V, AT, sz);
                atomnum = reshape(full(sum(V~=0, 2)), sz(1:end-1));
                mean_atomnum = full(mean(sum(V~=0, 2)));
                %set_beta_zero = true;

                if itn > 1

                    L1_old = full(sum(abs(V_old),2)) + reshape(sum(abs(abs(I) - abs(I_sparse_old)), 3), [prod(sz(1:end-1)) 1]);
                    L1_new = full(sum(abs(V),2)) + reshape(sum(abs(abs(I) - abs(I_sparse)), 3), [prod(sz(1:end-1)) 1]);

                    S_old = reshape(I_sparse_old, [prod(sz(1:end-1)), sz(end)]);
                    S = reshape(I_sparse, [prod(sz(1:end-1)), sz(end)]);

                    idx = L1_old < L1_new;

                    disp(['Not updating ' num2str(sum(L1_old < L1_new) / numel(L1_old)*100) '% of voxels'])

                    S(idx, :) = S_old(idx, :);
                    V(idx, :) = V_old(idx, :);

                    I_sparse = reshape(S, sz);

                end

                % scale up if L1 norm underestimated intensity
                I_sparse_vector = reshape(I_sparse, [prod(sz(1:end-1)), sz(end)]);
                I_vector = reshape(abs(I), [prod(sz(1:end-1)), sz(end)]);
                s = zeros(prod(sz(1:end-1)), 1);
                for i = 1:length(s)
                    s(i) = I_vector(i,:) / I_sparse_vector(i,:);
                end

                V = bsxfun(@times, V, s);
                I_sparse = inverse_basis_function(V, AT, sz);


            else
                % this is when L1 = 0
                I_sparse = zeros(size(I));
                V = sparse(prod(sz(1:end-1)), size(A,2));
            end




            % 2. get derivative of I
            % dCF/dI = d/dI { ||K - F(I)||^2 + lambda_2||I||_TV }

            [CF, I_grad] = CF_dictionary_2D(I, K_us, mask, samp_pdf, beta_TV, V, I_sparse, beta_L1, TE_mask, T2, AT);

            CF_begin = CF.CF;

            disp(['Cost function = ' num2str(CF.CF)])

            x_n = -I_grad;


            % 2. compute beta using Polak-Ribière formula
            %beta = (x_n(:)'*(x_n(:) - x_n_minus_one(:))) / (x_n_minus_one(:)'*x_n_minus_one(:) + eps);
            beta = (x_n(:)'*x_n(:)) / (x_n_minus_one(:)'*x_n_minus_one(:) + eps);


            beta = constrain(abs(beta), 0, inf); % provides automatic direction reset

            if set_beta_zero
                beta = 0;

            end

            beta_all(itn) = beta;


            disp(['Beta = ' num2str(beta_all(itn))])

            % 3. update the conjugate direction
            s_n = x_n + beta.*s_n_minus_one;


            % 4. perform a line search

            % only let cg alpha go between 1E-2 and 2
            cg_alpha = constrain(cg_alpha, 1E-6, 2);

            alpha = cg_alpha;

            min_line_its = 3;
            max_line_its = 20;


            %disp('Line iteration 0')
            I_potential = I + alpha*(s_n);% + V_s_n

            CF_potential = CF_dictionary_2D(I_potential, K_us, mask, samp_pdf, beta_TV, V, I_sparse, beta_L1, TE_mask, T2);


            % Do a line search
            alpha_line = [0 alpha];
            CF_line = [CF.CF CF_potential.CF];

            % for display purposes only
            Line_alpha = [0 alpha];
            Line_DC = [CF.DC CF_potential.DC];
            Line_L1 = [CF.L1 CF_potential.L1];
            Line_TV = [CF.TV CF_potential.TV];
            Line_CF = [CF.CF CF_potential.CF];


            best_CF = CF;
            best_alpha = 0;


            for line_its = 1:max_line_its

                %disp(['Line iteration ' num2str(line_its)])

                alpha = golden_ratio_search(alpha_line, CF_line);

                I_potential = I + alpha*(s_n ); %+ V_s_n

                % only call the first arguments - doesn't compute the residuals
                CF_potential = CF_dictionary_2D(I_potential, K_us, mask, samp_pdf, beta_TV, V, I_sparse, beta_L1, TE_mask, T2, AT);


                % for display purposes only
                Line_alpha = [Line_alpha alpha];
                Line_DC = [Line_DC CF_potential.DC];
                Line_L1 = [Line_L1 CF_potential.L1];
                Line_TV = [Line_TV CF_potential.TV];
                Line_CF = [Line_CF CF_potential.CF];

                figure(3);
                [val, idx] = sort(Line_alpha);
                subplot(2,2,1); plot(Line_alpha(idx), Line_DC(idx), '-*'); title('DC')
                subplot(2,2,2); plot(Line_alpha(idx), Line_L1(idx), '-*'); title('L1')
                subplot(2,2,3); plot(Line_alpha(idx), Line_TV(idx), '-*'); title('TV')
                subplot(2,2,4); scatter(Line_alpha(idx), Line_CF(idx), '*'); title('Line search cost')
                [fitres, min_alpha] = fit_polynomial(Line_alpha(idx), Line_CF(idx));
                hold on; plot(fitres, 'r'); hold off
                pause(0.1)

                % update the alpha and CF list
                alpha_line = [alpha_line alpha];
                CF_line = [CF_line CF_potential.CF];


                % if a new minimum has been found, update it
                if CF_potential.CF <= best_CF.CF;
                    best_CF = CF_potential;

                    best_alpha = alpha;
                end

                if (line_its >= min_line_its) && (best_alpha ~= 0) && (best_alpha ~= max(alpha_line))
                    break
                end
            end



            [fitres, min_alpha] = fit_polynomial(Line_alpha(idx), Line_CF(idx));

            if (min_alpha < 0) || (min_alpha > max(alpha_line))
                min_alpha = best_alpha;
            end

            if min_alpha == 0
                min_alpha = 1E-3;
            end

            cg_alpha = min_alpha;
            disp(['Minimum found with alpha = ' num2str(cg_alpha)])

            % update I
            I = I + min_alpha*(s_n);% + V_s_n);


    %         if rem(itn, 20) == 0
    %             set_beta_zero = true;
    %         else
    %             set_beta_zero = false;
    %         end


            % handle the previous iteration
            CF_old = CF;
            x_n_minus_one = x_n;
            s_n_minus_one = s_n;
            %V_x_n_minus_one = V_x_n;
            %V_s_n_minus_one = V_s_n;


            figure(1); imagesc_nD(abs(I)); drawnow

            k = fft2c(I);

            figure(2); imagesc(log(abs(k(:,:,1))), [-6 6]); drawnow

            figure(4); plot(squeeze(abs(I(:,end/2,5)))); drawnow


            [ DT, FA, ADC, VectorField, EigVals] = fit_DT( abs(I)+eps, bval);
            HA = compute_helix_angle_on_image(rad_vec, circ_vec, long_vec, VectorField(:,:,[2 3 1],1));

            HA(HA > (HA_orig + 45)) = HA(HA > (HA_orig + 45)) - 180;
            HA(HA < (HA_orig - 45)) = HA(HA < (HA_orig - 45)) + 180;


            T = Tissue;%imerode(Tissue, ones(3));

            FA_rmse(itn) = sqrt(mean((FA_orig(T) - FA(T)).^2));
            ADC_rmse(itn) = sqrt(mean((ADC_orig(T) - ADC(T)).^2));
            HA_rmse(itn) = sqrt(mean((HA_orig(T) - HA(T)).^2));
            
            FA_mean(itn) = mean(FA(T));
            ADC_mean(itn) = mean(ADC(T));
            FA_std(itn) = std(FA(T));
            ADC_std(itn) = std(ADC(T));

            figure(6); subplot(1,3,1); plot(FA_rmse); title('FA RMSE'); % 0.0243
            subplot(1,3,2); plot(ADC_rmse * 1000); title('ADC RMSE'); % 0.0502
            subplot(1,3,3); plot(HA_rmse); title('HA RMSE'); % 3.138
            drawnow

            figure(7); imagesc((HA - HA_orig).*T);

    %         if best_CF.DC < (2.7E4 / accel_factor) % DC too small, overfitting noise
    %             beta_TV = beta_TV * 1.1
    %             beta_L1 = beta_L1 * 1.1 
    %         else
    %             beta_TV = beta_TV / 1.1
    %             beta_L1 = beta_L1 / 1.1
    %         end

            DC_all(itn) = best_CF.DC;
            TV_all(itn) = best_CF.TV;
            L1_all(itn) = best_CF.L1;
            figure(8); plot(DC_all, 'b');
            hold on; plot([1 itn], [2.7E4 / accel_factor, 2.7E4 / accel_factor], '--r')
            hold off
        end

        FA_rmse_accel(reps, accel_factor-1) = FA_rmse(itn);
        ADC_rmse_accel(reps, accel_factor-1) = ADC_rmse(itn);
        HA_rmse_accel(reps, accel_factor-1) = HA_rmse(itn);
        FA_mean_accel(reps, accel_factor-1) = FA_mean(itn);
        ADC_mean_accel(reps, accel_factor-1) = ADC_mean(itn);
        FA_std_accel(reps, accel_factor-1) = FA_std(itn);
        ADC_std_accel(reps, accel_factor-1) = ADC_std(itn);
        
        
    end
end




retro_FA = mean(FA_rmse_accel);
retro_ADC = mean(ADC_rmse_accel);
retro_HA = sort(mean(HA_rmse_accel));


prosp_FA = mean(FA_rmse_accel);
prosp_ADC = mean(ADC_rmse_accel);
prosp_HA = sort(mean(HA_rmse_accel));






