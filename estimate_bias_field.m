function bias = estimate_bias_field(im, ridx, cidx, sidx, rad)
%
% As part of Becky Burton's paper for Euroecho 2014, I worked on this
% algorithm to estimate the bias field by splitting the image into blocks,
% and then in each block smoothing the histogram until we get 1 peak
% (background only) or two peaks (background and foreground object). The
% background is estimated as the mode of the background peak. But this
% function is not quite finished. The histogram smoothing can get trapped
% into an infinite loop where the number of peaks doesn't decrease. A
% similar idea I implemented more properly in FiltersToolbox/amrf_seg.m,
% but that one uses the background and foreground peaks to do Markov Random
% Field segmentation.

% initialize array to contain the estimated background intensity
bias_sample = nan(length(ridx), length(cidx), length(sidx));

% create the image blocks that are going to be used to estimate each bias
% value
bias_block = cell(length(ridx), length(cidx), length(sidx));
tic
for R = 1:length(ridx)
    for C = 1:length(cidx)
        for S = 1:length(sidx)

            bias_block{R, C, S} = im(...
                max(1, ridx(R)-rad(1)):min(size(im, 1), ridx(R)+rad(1)), ...
                max(1, cidx(C)-rad(2)):min(size(im, 2), cidx(C)+rad(2)), ...
                max(1, sidx(S)-rad(3)):min(size(im, 3), sidx(S)+rad(3)) ...
                );
            sz = size(bias_block{R, C, S});
            bias_block{R, C, S} ...
                = bias_block{R, C, S}(bias_block{R, C, S} > 0);
            
%             % avoid wasting time estimating the bias field for sampling
%             % points that are mostly in a masked-out region
%             if (numel(bias_block{R, C, S}) < prod(sz)*.75)
%                 bias_block{R, C, S} = [];
%             end
   
        toc
        end
    end
end
disp('Blocks created')

% process the blocks in parallel, computing the bias value in each block
tic
for R = 1:length(ridx)
    for C = 1:length(cidx)
        for S = 1:length(sidx)
            [R, C, S]
            bias_sample(R, C, S) = estimate_local_bias(...
                bias_block{R, C, S});
            toc
        end
    end
end
disp('Bias estimated')

% upsample the bias field values to obtain a bias field for the original
% image. To avoid running out of memory, we are going to take each
% bias_sample and its 8 neighbours, and use that to interpolate the part of
% the image they comprise
bias = zeros(size(im));
for R = 2:length(ridx)-1
    for C = 2:length(cidx)-1
        for S = 2:length(sidx)-1
            
            % coarse grid with bias_sample and its 8 neighbours
            [gr, gc, gs] ...
                = ndgrid(ridx(R-1:R+1), cidx(C-1:C+1), sidx(S-1:S+1));
            gbias_sample = bias_sample(R-1:R+1, C-1:C+1, S-1:S+1);
            
            % fine grid of the original image
            [ggr, ggc, ggs] = ndgrid(...
                ridx(R-1):ridx(R+1), ...
                cidx(C-1):cidx(C+1), ...
                sidx(S-1):sidx(S+1) ...
                );
            
            % thin-plate spline interpolation of the fine grid from the
            % coarser grid (ignoring bias field estimates that are NaNs)
            idx = ~isnan(gbias_sample(:));
            ggbias_sample = pts_tps_map([gr(idx), gc(idx), gs(idx)], ...
                gbias_sample(idx), [ggr(:), ggc(:), ggs(:)]);
            ggbias_sample = reshape(ggbias_sample, ...
                ridx(R+1) - ridx(R-1) + 1, ...
                cidx(C+1) - cidx(C-1) + 1, ...
                sidx(S+1) - sidx(S-1) + 1 ...
                );
            
            % save interpolation result into output array
            bias(...
                ridx(R-1):ridx(R+1), ...
                cidx(C-1):cidx(C+1), ...
                sidx(S-1):sidx(S+1) ...
                ) = ggbias_sample;
            
%             % compute interpolant
%             biasf = TriScatteredInterp(gr(:), gc(:), gs(:), gbias_sample(:), 'linear');
% 
%             % interpolate bias field values
%             tic
%             aux =  biasf(gr(:), gc(:), gs(:));
%             toc
%             bias(...
%                 ridx(R-1):ridx(R+1), ...
%                 cidx(C-1):cidx(C+1), ...
%                 sidx(S-1):sidx(S+1) ...
%                 ) = biasf(gr(:), gc(:), gs(:));
            
        end
    end
end





end


function mubak = estimate_local_bias(im)

if (isempty(im))
    mubak = nan;
    return
end

% compute histogram of the intensity values
[fhist, xhist] = hist(im, 500);
fhist = fhist / sum(fhist);

% smooth the histogram until we find just 1 or 2 peaks
npks = Inf;
tol = 0.5e-10;
while (1)
    
    % increase the smoothing parameter
    tol = tol * 2;
    [~, fhist2] = spaps(xhist, fhist, tol);
    
    % find the peaks
    [pks, loc] = findpeaks(fhist2, 'minpeakheight', 0.5e-3);
    
    % number of peaks found
    npks = length(pks);
    
    % stop conditions
    if ((npks == 1) ...
            || ((npks == 2) && (abs(diff(loc)) >= 50)))
        break
    end
    
end

% deal with the number of peaks
if (npks == 2)
    
    % we have background and tissue. The background is lighter. Sort the
    % peaks
    [loc, idx] = sort(loc, 'descend');
    pks = pks(idx);
    
    % DEBUG
    hold off
    plot(xhist, fhist, 'b')
    hold on
    plot(xhist, fhist2, 'r')
    plot(xhist(loc(1))*[1 1], [0 max(fhist)], 'r')
    plot(xhist(loc(2))*[1 1], [0 max(fhist)], 'r')

    mubak = xhist(max(loc));
    
elseif (npks == 1)
    
    % we assume there's only background (although note that this could be a
    % case with only tissue)
    mubak = xhist(loc);
    
%     % DEBUG
%     hold off
%     plot(xhist, fhist, 'b')
%     hold on
%     plot(xhist, fhist2, 'r')
%     plot(xhist(loc)*[1 1], [0 max(fhist)], 'r')
    
else

    error('Assertion fail: The histogram has no peaks')
    
end

% % compute mixture of Gaussians models for 1 Gaussian and 2 Gaussians
% warning('off', 'stats:gmdistribution:FailedToConverge')
% obj1 = gmdistribution.fit(double(im), 1);
% obj2 = gmdistribution.fit(double(im), 2);
% warning('on', 'stats:gmdistribution:FailedToConverge')
% 
% % select which model fits the data better
% if (obj1.Converged && obj2.Converged)
%     
%     % if both Gaussian mixture models (1 or 2 Gaussians) are feasible,
%     % the best one is the one with the lowest Bayesian Information
%     % Criterion (BIC)
%     if (obj1.BIC <= obj2.BIC)
%         nobj = 1;
%     else
%         nobj = 2;
%     end
%     
% elseif (obj1.Converged)
%     
%     % if only the 1-Gaussian model converges, that's the one we select
%     nobj = 1;
%     
% elseif (obj2.Converged)
%     
%     % if only the 2-Gaussian model converges, that's the one we select
%     nobj = 2;
%     
% else
%     
%     % if neither the 1 or 2-Gaussian models converge, throw an error
%     error('The image cannot be modelled with 1 or 2 Gaussian pdfs')
%     
% end
% 
% % even if the 2-Gaussian model is a good fit to the data, in practice
% % maybe what we have is only tissue or background, but with enough
% % variability to trick us into thinking 2 Gaussians are a better fit.
% % If we cannot find a good separation between the 2 Gaussians, we
% % assume that the 1 Gaussian model is better
% if (nobj == 2)
%     
%     
%     % get mixture of Gaussians parameters
%     [mutis, idx] = min(obj2.mu);
%     vartis = obj2.Sigma(idx);
%     [mubak, idx] = max(obj2.mu);
%     varbak = obj2.Sigma(idx);
%     
%     % DEBUG: create Gaussian curves for display purposes
%     %     hold off
%     %     [fhist, xout] = hist(im, 500);
%     %     ftis = normpdf(xout, mutis, sqrt(vartis));
%     %     fbak = normpdf(xout, mubak, sqrt(varbak));
%     %     hold off
%     %     plot(xout, fhist/sum(fhist), 'b')
%     %     hold on
%     %     plot(xout, ftis/sum(ftis), 'r')
%     %     plot(xout, fbak/sum(fbak), 'g')
%     
%     % compute intersection points between two gaussians
%     thr = intersect_gaussians(mutis, mubak, sqrt(vartis), sqrt(varbak));
%     
%     % keep the one that is between both mean values
%     thr = thr(thr > mutis & thr < mubak);
%     
%     if (isempty(thr) || isnan(thr))
%         % if there's no threshold between the two Gaussians, the quality
%         % measure is NaN
%         q = nan;
%         
%         % we cannot get a separation between the two Gaussians, so we
%         % assume that intensities actually follow a single
%         nobj = 1;
%     else
%         % quality parameter is a measure of gaussian overlap
%         q = 1 - normcdf(2*mutis-thr, mutis, sqrt(vartis))...
%             - normcdf(thr, mubak, sqrt(varbak));
%     end
%     
% end
% 
% 
% % in a 1-Gaussian model, we are going to assume that the block contains
% % only background voxels
% if (nobj == 1)
%     
%     mubak = obj1.mu;
%     
% end

end
