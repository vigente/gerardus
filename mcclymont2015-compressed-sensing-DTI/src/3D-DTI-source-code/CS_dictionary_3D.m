% this version works with Irivn's FSE data

close all
clear all
clc

scan_dirs = {'REDACTED'};


for sample_counter = 1:5
    
    data_dir = scan_dirs{sample_counter};
    

    file_dir = {'REDACTED'};%6x

    ground_truth_dir = 'REDACTED';

    noise_dir = 'REDACTED';

    noise = getfield(load_untouch_nii([data_dir, noise_dir, 'image_mag.nii']), 'img') ...
        .* exp(1i * getfield(load_untouch_nii([data_dir, noise_dir, 'image_phase.nii']), 'img'));

    noise_gain_dB = getfield(search_procpar([data_dir, noise_dir, 'procpar'], 'gain'), 'gain');

    % load the prescan power array
    load([data_dir, 'REDACTED'])
    db = gain.gain;
    scale = gain.scale; 
    
    for retro = 0:1

        for file_counter = 1:length(file_dir)

            % load the undersampled k-space
            [k_space_undersampled, TE_mask] = sort_cs([data_dir, file_dir{file_counter}]);
            sz = size(k_space_undersampled);

            % load TE
            pp = search_procpar('REDACTED', 'te', 'esp');
            TE = pp.te + (0:7) * pp.esp;


            if retro

                % retrospective under-sampling

                TE_mask(:) = 0;

                for i = 1:length(TE)

                    inds = floor((0:159) / 10);
                    tops = inds == (8 - i);
                    bottoms = inds == (7 + i);

                    TE_mask(:,tops,:,:) = TE(i);
                    TE_mask(:,bottoms,:,:) = TE(i);

                end

            else

                % prospective undersampling
                for i = 1:8
                    TE_mask(TE_mask==i) = TE(i);
                end
            end


            k_space_undersampled = k_space_undersampled(end:-1:1, end:-1:1, end:-1:1, :) * sqrt(prod(sz(1:3)));

            % get the b values
            bval = procpar_to_bval([data_dir, file_dir{file_counter}, 'procpar']);

            % re-order the sequence
            dwi_order = [1, 9, 17, 25, setdiff(1:34, [1, 9, 17, 25])];
            k_space_undersampled = k_space_undersampled(:,:,:,dwi_order);
            bval = bval(:,:,dwi_order);

            % load the fully sampled image
            ground_truth = getfield(load_untouch_nii([data_dir, ground_truth_dir, 'image_mag.nii']), 'img');
            dwi_order_gt = [1, 9, 17, 25, setdiff(1:35, [1, 9, 17, 25, 33])];
            ground_truth = ground_truth(:,:,:,dwi_order_gt);

            % get the mask of sampled points
            mask = k_space_undersampled ~= 0;

            if retro
                % test
                %mask(:) = 1;
                % 
                k_space_undersampled = fft3c(ground_truth) .* mask;
            end

            % get the data (zero filled)
            I_orig = abs(ifft3c(k_space_undersampled));

            % get a rough segmentation
            M = (mean(ground_truth(:,:,:,1:4),4) > 10000) & (mean(ground_truth(:,:,:,1:4),4) ./ mean(ground_truth(:,:,:,5:end),4) < 6);

            % normalise the signal so that the noise is equal for all scans
            DTI_gain_dB = getfield(search_procpar([data_dir, file_dir{file_counter}, 'procpar'], 'gain'), 'gain');
            DTI_gain_dB = DTI_gain_dB(dwi_order);

            noise_level = zeros(1, size(noise,4));
            for i = 1:size(noise,4)
                nn = real(noise(:,:,:,i));
                noise_level(i) = std(nn(:));
            end

            noise_normaliser = zeros(1,1,1,size(I_orig,4));
            noise_normaliser(1:4) = noise_level(1);
            noise_normaliser(5:end) = noise_level(2);

            k_space_undersampled = bsxfun(@rdivide, k_space_undersampled, noise_normaliser); 



        %     % load the dictionary - trained on older data
        %     load('REDACTED')
        %     AT = A';
        %     
        %     % this is temporary, until I have properly trained dictionaries
        %     AT = bsxfun(@rdivide, AT, squeeze(noise_normaliser)');
        %     AT = bsxfun(@rdivide, AT, sqrt(sum(AT.^2, 2)));
        %     A = AT';
        %     
            
            load(['REDACTED' num2str(sample_counter) '.mat'])

            % prospective undersampling, temperature should be stable by now
            if ~retro        
                A(1:4,:) = repmat(mean(A(2:4,:),1), [4 1]);
            end

            AT = A';


            % load the T2 array
            cd([data_dir, 'REDACTED'])
            Imag = load_untouch_nii('image_mag.nii');
            Imag = Imag.img;

            Iphase = load_untouch_nii('image_phase.nii');
            Iphase = Iphase.img;

            ME3D = Imag .* exp(1i * Iphase);

            sz = size(ME3D);

            % fit T2
            I_vector = reshape(abs(double(ME3D)), [prod(sz(1:3)), sz(4)]);

            TE_ME3D = 1/1000*getfield(search_procpar('procpar', 'TE'), 'TE');

            TE_mat = [ones(size(TE_ME3D)); TE_ME3D]';

            % linear fit
            %AB = pinv(TE_mat) * log(I_vector)'; % same as (log(I_vector') / TE_mat')'

            % choose a representative signal from linear fit (anything but noise)
            %idx = sub2ind(sz(1:3), 100, 55, 80);
            %d = AB(:, idx);
            d = [10.9179; -40.6271];

            % weighting matrix (N x N)
            W = eye(length(TE_mat));
            W(W == 1) = exp(TE_mat * d);

            % weighted fit
            AB = pinv(TE_mat' * (W.^2) * TE_mat) * TE_mat' *(W.^2) * log(I_vector)';

            %A = exp(reshape(AB(1,:), sz(1:3)));
            T2 = -1 ./ reshape(AB(2,:), sz(1:3));

            T2(abs(ME3D(:,:,:,1)) < 5000) = inf;


        %     sz = size(I_orig);
        %     
        %     S = reshape(I_orig, [prod(sz(1:end-1)), sz(end)]);
        %     S = S(M(:), :);
        %     
        %     
        %     S(:, 1:4) = repmat(mean(S(:,1:4), 2), [1, 4]);
        %     
        %     
        %     % generate the starting conditions for the dictionary
        %     D = generate_tensor_dictionary(0.6, 1.3E-3, 258, bval)';
        %     D_normaliser = sqrt(sum(D.^2, 1));
        %     D_norm = D ./ repmat(D_normaliser, [size(D,1), 1]);
        %     
        %     
        %     param.iter = 200;
        %     param.K = 100; % 258
        %     param.lambda=0.58^2 * sqrt(34);
        %     param.mode=1;
        %     param.posAlpha=1;
        %     param.posD=1;
        %     param.clean = 1;
        %     param.batchsize = 512;
        %     param.verbose = 1;
        % 
        %     % train the model
        %     [U, model] = mexTrainDL(S',param);
        % 
        %     A = U;
        % 
        %     % ensure dictionary is L2 normalised
        %     l_norm = 2;
        %     A = A ./ repmat((sum(A.^l_norm ,1).^(1/l_norm)), [size(A,1), 1]);
        %     AT = A';
        % 
        %     return

            sz = size(I_orig);

            samp_pdf = zeros(sz);

            accel_factor = round(numel(k_space_undersampled) / sum(k_space_undersampled(:) ~= 0));

            central_region = 0.15;

            num_slices = sz(2);
            num_PE = sz(3);
            num_echos = 8;

            for i = 1:size(mask,4)

                [mask_2D, PE1, PE2, samp_pdf_dir] = generate_FSE_DTI_mask(accel_factor, central_region, num_slices, num_PE, num_echos);
                samp_pdf(:,:,:,i) = repmat(permute(samp_pdf_dir, [3 1 2]), [sz(1) 1 1]);

            end

            disp('Running sliding window')
            k_zero = DTI_sliding_window(k_space_undersampled(:,:,:,1:4), bval(:,:,1:4));

            k_high = DTI_sliding_window(k_space_undersampled(:,:,:,5:end), bval(:,:,5:end));

            % try and get any other points, even if they are from bhigh
            k_zero = DTI_sliding_window(cat(4, k_zero, k_high*2), bval);
            k_zero = k_zero(:,:,:,1:4);
            
            %I = bsxfun(@rdivide, ifft3c(cat(4, k_zero, k_high)), exp(- TE(1) / T2));
            I = ifft3c(cat(4, k_zero, k_high));


            clear k_high
            clear k_zero


            I_phase = angle(I);

            [ DT_orig, FA_orig, ADC_orig, VectorF_orig ] = fit_DT( abs(ground_truth(100, :, :,2:end)), ...
                                                                   bval(:,:,2:end), 0 );

            FA_orig = squeeze(FA_orig);
            ADC_orig = squeeze(ADC_orig);


            % In retrospective recons: TV = 1E-3, L1 = 1;

            % Mani paper: b1 = 1e-4 or 1e-5, b2 = 1e-5
            beta_TV = 1E-4;
            % sample 1 was done with 1E-3 and 1
            beta_L1 = 0.1;   % L1


            % test for fully sampled image
            %beta_TV = 1E-3;
            %beta_L1 = 0;
            % delete this line
        % for ground truth T2w compensation, TV = 1E-3 was still not enough.


            x_n_minus_one = zeros(size(I));
            s_n_minus_one = zeros(size(I));
            %V_x_n_minus_one = zeros(size(I));
            %V_s_n_minus_one = zeros(size(I));

            cg_alpha = 0.5;
            CF_old.CF = inf;

            disp_ind = sub2ind(sz(1:3), 100, 90, 60);

            set_beta_zero = false;

            sparsity_error = 1;

            for itn = 1:70

                disp(['Iteration ' num2str(itn)])

                % 1. get best estimate of V
                % V* = argmin_V { ||I - D*V||_2^2 + lambda_1||V||_1 }

                figure(31); subplot(321);hold off; plot(squeeze(abs(I(100,90, 60, :))));

                if itn > 1
        %                 sparsity_error = min(sparsity_error, ...
        %                     mean(mean(mean(sum((abs(I) - abs(I_sparse)).^2, 4)))) / sz(end));
        % ^ this worked, but I want to see if I can do it without the min
                    sparsity_error = mean(mean(mean(sum((abs(I) - abs(I_sparse)).^2, 4)))) / sz(end);
                end


                if beta_L1 > 0

                    if itn > 1
                        V_old = V;
                        I_sparse_old = I_sparse;
                    end

                    disp(['Running forward transform, allowable error = ' num2str(sparsity_error)])
                    V = forward_basis_function(abs(I), AT, sz, sparsity_error);

                    I_sparse = inverse_basis_function(V, AT, sz);
                    atomnum = reshape(full(sum(V~=0, 2)), sz(1:end-1));
                    mean_atomnum = full(mean(sum(V~=0, 2)));
                    %set_beta_zero = true;

                    if itn > 1

                        L1_old = full(sum(abs(V_old),2)) + reshape(sum(abs(abs(I) - abs(I_sparse_old)), 4), [prod(sz(1:3)) 1]);
                        L1_new = full(sum(abs(V),2)) + reshape(sum(abs(abs(I) - abs(I_sparse)), 4), [prod(sz(1:3)) 1]);

                        S_old = reshape(I_sparse_old, [prod(sz(1:3)), sz(4)]);
                        S = reshape(I_sparse, [prod(sz(1:3)), sz(4)]);

                        idx = L1_old < L1_new;

                        disp(['Not updating ' num2str(sum(L1_old < L1_new) / numel(L1_old)*100) '% of voxels'])

                        S(idx, :) = S_old(idx, :);
                        V(idx, :) = V_old(idx, :);

                        I_sparse = reshape(S, sz);

                    end

                    figure(31); subplot(322); plot(V(disp_ind, :));


                    figure(40); imagesc(squeeze(atomnum(100, :, :))); title(['Mean atomnum = ' num2str(mean_atomnum) ]); colorbar

                    figure(31); subplot(321); hold on; plot(squeeze(I_sparse(100, 90, 60, :)), 'r'); 
                    plot(squeeze(abs(ground_truth(100, 90, 60,:))) ./ noise_normaliser(:), 'g')
                    pause(0.1)

                else
                    % this is when L1 = 0
                    I_sparse = zeros(size(I));
                    V = sparse(prod(sz(1:3)), sz(4));
                end



                %V_x_n = - (abs(I) - I_sparse).*exp(1i*I_phase); % "derivative" of V

                %figure(31); subplot(323); plot(squeeze(real(V_x_n(90, 60,:) .* exp(-1i*I_phase(90,60,:))))); title('V_x_n magnitude')


                %beta_V =  (V_x_n(:)'*(V_x_n(:))) / (V_x_n_minus_one(:)'*V_x_n_minus_one(:) + eps);

                %beta_V = constrain(abs(beta_V), 0, 0.95);



                %disp(['V beta = ' num2str(beta_V)])


                % 2. get derivative of I
                % dCF/dI = d/dI { ||K - F(I)||^2 + lambda_2||I||_TV }

                [CF, I_grad] = CF_dictionary_3D(I, k_space_undersampled, mask, samp_pdf, beta_TV, V, I_sparse, beta_L1, TE_mask, T2, AT);

                CF_begin = CF.CF;

                disp(['Cost function = ' num2str(CF.CF)])

                x_n = -I_grad;

                figure(31); subplot(324); plot(squeeze(real(x_n(100,90, 60, :).*exp(-1i*I_phase(100,90,60,:))))); title('x_n')

                % 2. compute beta using Polak-Ribière formula
                %beta = (x_n(:)'*(x_n(:) - x_n_minus_one(:))) / (x_n_minus_one(:)'*x_n_minus_one(:) + eps);
                beta = (x_n(:)'*x_n(:)) / (x_n_minus_one(:)'*x_n_minus_one(:) + eps);

                % I think it neeeds to be this one, because the scaling of the
                % TV is not certain, and it needs time to build up. If we use
                % the other one, the x_n and x_n_minus_one are quite similar,
                % so we get a low beta, and the gradient can't accumulate.


                beta = constrain(abs(beta), 0, inf); % provides automatic direction reset

                if set_beta_zero
                    beta = 0;
                    cg_alpha = 0.5;
                end

                beta_all(itn) = beta;


                disp(['Beta = ' num2str(beta_all(itn))])

                % 3. update the conjugate direction
                s_n = x_n + beta.*s_n_minus_one;

                figure(31); subplot(325); hold off; plot(squeeze(real(s_n(100,90, 60, :).*exp(-1i*I_phase(100,90,60,:))))); title('s_n')


                %V_s_n = V_x_n + beta* V_s_n_minus_one;

                %figure(31); subplot(325); hold on; plot(squeeze(real(V_s_n(90, 60, :) .*exp(-1i*I_phase(90,60,:)))), 'r');

                %figure(31); subplot(326); plot(squeeze(real((V_s_n(90, 60, :) + s_n(90, 60, :)) .*exp(-1i*I_phase(90,60,:)))), 'r'); title('S_n sum')


                % 4. perform a line search

                % only let cg alpha go between 1E-2 and 2
                cg_alpha = constrain(cg_alpha, 0, 2);

                alpha = cg_alpha;

                min_line_its = 3;
                max_line_its = 10;


                %disp('Line iteration 0')
                I_potential = I + alpha*(s_n);% + V_s_n

                CF_potential = CF_dictionary_3D(I_potential, k_space_undersampled, mask, samp_pdf, beta_TV, V, I_sparse, beta_L1, TE_mask, T2);


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
                    CF_potential = CF_dictionary_3D(I_potential, k_space_undersampled, mask, samp_pdf, beta_TV, V, I_sparse, beta_L1, TE_mask, T2);


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

                I_phase = angle(I);



                % handle the previous iteration
                CF_old = CF;
                x_n_minus_one = x_n;
                s_n_minus_one = s_n;
                %V_x_n_minus_one = V_x_n;
                %V_s_n_minus_one = V_s_n;


                k_space_soln = fft3c(I);




                figure(5)
                subplot(4,2,[1 3]); imagesc(squeeze(log(abs(k_space_undersampled(100, :, :, 2)))), [-2 8]); title('Undersampled K-space')
                subplot(4,2,[2 4]); imagesc(squeeze(log(abs(k_space_soln(100, :, :, 2)))), [-2 8]); title('Reconstructed K-space')
                DC_all(itn) = CF.DC;
                TV_all(itn) = CF.TV;
                CF_all(itn) = CF.CF;
                L1_all(itn) = CF.L1;
                subplot(4,3, 7:8); plot(CF_all,'r');
                subplot(4,3,9); plot(1, 'r'); hold on; plot(1,'b'); plot(1,'g'); plot(1, 'k'); legend('CF', ' DC', 'TV', 'L1')
                subplot(4,3,10); plot(DC_all, 'b');
                subplot(4,3,11); plot(TV_all, 'g');
                subplot(4,3,12); plot(L1_all, 'k')


                I_denorm = bsxfun(@times, I, noise_normaliser);

                % simulate the k-space of the ground truth
                k_space = zeros(size(k_space_undersampled));
                for i = 1:length(TE)

                    inds = floor((0:159) / 10);
                    tops = inds == (8 - i);
                    bottoms = inds == (7 + i);

                    % get an idea of what the decayed image should be
                    I_decayed = bsxfun(@times, I_denorm, exp(- (TE(i)-TE(1)) / T2));

                    % get the fourier transform of the decayed image (a portion of it
                    % will be used later)
                    K_decayed = fft3c(I_decayed);

                    % put the right values into k-space
                    k_space(:,tops,:,:) = K_decayed(:,tops,:,:);
                    k_space(:,bottoms,:,:) = K_decayed(:,bottoms,:,:);


                end

                I_recon = ifft3c(k_space);


                [ DT, FA, ADC, VectorF ] = fit_DT( abs(I_recon(100, :, :,2:end)), bval(:,:,2:end), 0 );

                figure(30); 
                subplot(451); imagesc(squeeze(FA_orig), [0 0.5]);  colormap gray;
                subplot(452); imagesc(squeeze(FA), [0 0.5]);  colormap gray;
                subplot(453); imagesc(squeeze(FA) - FA_orig, [-0.1 0.1]);  colormap jet;

                subplot(456); imagesc(squeeze(ADC_orig), [0 3E-3]); colormap gray;
                subplot(457); imagesc(squeeze(ADC), [0 3E-3]); colormap gray;
                subplot(458); imagesc(squeeze(ADC) - ADC_orig, [-0.3E-3 0.3E-3]); colormap jet

                subplot(4,5,11); imagesc(squeeze(abs(ground_truth(100,:,:,2))), [0 6E4]); colormap gray;
                subplot(4,5,12); imagesc(squeeze(abs(I_recon(100,:,:,2))), [0 6E4]); colormap gray;
                subplot(4,5,13); imagesc(squeeze(abs(I_recon(100,:,:,2)) - abs(ground_truth(100,:,:,2))), [-5E3 5E3]); colormap jet

                subplot(4,5,16); imagesc(squeeze(abs(ground_truth(100,:,:,10))), [0 3E4]); colormap gray;
                subplot(4,5,17); imagesc(squeeze(abs(I_recon(100,:,:,10))), [0 3E4]); colormap gray;
                subplot(4,5,18); imagesc(squeeze(abs(I_recon(100,:,:,10)) - abs(ground_truth(100,:,:,10))), [-3E3 3E3]); colormap jet

                pause(0.1)

                FA = squeeze(FA);
                ADC = squeeze(ADC);

                DIFF = FA - FA_orig;

                M_slice = M(100, :, :);

                NRMSE_FA(itn) = sqrt(mean(DIFF(M_slice).^2));


                DIFF = ADC - ADC_orig;

                NRMSE_ADC(itn) = sqrt(mean(DIFF(M_slice).^2));

                DIFF = squeeze(abs(ground_truth(100,:,:,10)) - abs(I_recon(100,:,:,10)));

                NRMSE_I(itn) = sqrt(mean(DIFF(M_slice).^2));


                [ DT, FA, ADC, VectorF ] = fit_DT( abs(I_denorm(99, :, :,2:end)), bval(:,:,2:end), 0 );

                subplot(454); imagesc(squeeze(FA), [0 0.5]); 
                subplot(459); imagesc(squeeze(ADC), [0 3E-3]); 
                subplot(4,5,14); imagesc(squeeze(abs(I_denorm(100,:,:,2))), [0 6E4]);
                subplot(4,5,19); imagesc(squeeze(abs(I_denorm(100,:,:,10))), [0 3E4]);

                subplot(455); imagesc(squeeze(FA) - FA_orig, [-0.1 0.1]);  colormap jet;
                subplot(4,5,10); imagesc(squeeze(ADC) - ADC_orig, [-0.3E-3 0.3E-3]); colormap jet
                subplot(4,5,15); imagesc(squeeze(abs(I_denorm(100,:,:,2)) - abs(ground_truth(100,:,:,2))), [-5E3 5E3]); colormap jet
                subplot(4,5,20); imagesc(squeeze(abs(I_denorm(100,:,:,10)) - abs(ground_truth(100,:,:,10))), [-3E3 3E3]); colormap jet



                figure(38); subplot(231); plot(NRMSE_FA); title('RMSE FA')
                subplot(232); plot(NRMSE_ADC); title('RMSE ADC')
                subplot(233); plot(NRMSE_I); title('RMSE I')

                beta_all(1) = 0;
                subplot(2,3,4); plot(beta_all); title('Beta')

        %                 subplot(234); scatter(FA_orig(M), FA(M), '.'); hold on; plot([0.05 0.4], [0.05 0.4]); hold off;
        %                 
        %                 subplot(234); scatter(ADC_orig(M), (M), '.'); hold on; plot([0.05 0.4], [0.05 0.4]); hold off;
        %                
        %                 subplot(234); scatter(FA_orig(M), FA(M), '.'); hold on; plot([0.05 0.4], [0.05 0.4]); hold off;



            end


            I = bsxfun(@times, I, noise_normaliser);

           
            if retro
                rr = '_retro_';
            else
                rr = '_prosp_';
            end


            save_name = ['Sample_' num2str(sample_counter) ...
                            rr ...
                        'acc_' num2str(file_counter+1) ...
                           '.mat'];

            disp(['Saving ' save_name])
                 
            % save the file
            cd('REDACTED')
       
            save(save_name, 'I')
            
            pause(10);
        end
    end
end