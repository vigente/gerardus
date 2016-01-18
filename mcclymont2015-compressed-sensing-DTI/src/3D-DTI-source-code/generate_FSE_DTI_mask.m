function [mask, PE1, PE2, P] = generate_FSE_DTI_mask(accel_factor, central_region, num_slices, num_PE, num_echos)

% close all
% clear all
% clc
% 
% 
% accel_factor = 10;
% central_region = 0.1;
% 
% num_slices = 192;%160
% 
% num_PE = 192;
% 
% num_echos = 8;

poly_order = accel_factor+1;

P = genPDF([num_slices num_PE], poly_order, 1/accel_factor, 2, central_region, 0);

% Get the pdf along the central line of k-space
P_line = P(end/2, 1:end/2);
% cumulative density function - we want this to be even between the four
% bands
P_cum = cumsum(P_line) / sum(P_line);
P_cum = P_cum - min(P_cum);


% define the bands
bands = num_echos + 1 - round(P_cum*num_echos + 0.5);
bands = constrain(bands, 1, num_echos);

% it is possible that a band will be completely skipped out - check for
% this
for i = num_echos:-1:1
    b = bands == i;
    if sum(b) == 0
        % find the index of the one below it
        idx = find(bands == i+1, 1, 'last');
        bands(1:idx-1) = bands(2:idx); % shift them all
        bands(idx) = i; % add in the missing index
    end
end

% it is also possible that a higher band will have a smaller band width
% than a lower band because of rounding - we don't want this either


for i = 2:num_echos
    b = bands == i;
    b_lower = bands == i-1;
    while sum(b) < sum(b_lower)
        
        idx = find(bands > i, 1, 'last');
        bands(idx) = i; % shift them all
        
        b = bands == i;
        b_higher = bands == i+1;
    end
end
        


min_band_width = sum(bands == 1);

% how many shots does each slice get?
slice_P = sum(P,2);
fully_sampled_region = zeros(size(slice_P));
min_idx = round(length(slice_P)/2 - central_region/2*length(slice_P));
max_idx = round(length(slice_P)/2 + central_region/2*length(slice_P));
fully_sampled_region(min_idx:max_idx) = 1;
shots_per_slice = pdf_to_mask(slice_P/num_echos/2, ones(size(slice_P))*min_band_width, fully_sampled_region); % /2 for half of k-space
shots_per_slice_top = shots_per_slice;

% given n shots, how to sample band 1, 2 etc?
N_top = zeros([num_PE/2, num_slices]);
for i = 1:num_slices
    
    for b = 1:num_echos
        % band b
        norm_prob = P_line .* (bands == b);
        norm_prob = norm_prob / sum(norm_prob);
        
        N = pdf_to_mask(norm_prob * shots_per_slice(i), ones(size(norm_prob)));
        
        N_top(:,i) = N_top(:,i) + N';
    end
   
end

shots_per_slice = pdf_to_mask(slice_P/num_echos/2, ones(size(slice_P))*min_band_width, fully_sampled_region); % /2 for half of k-space
% shots_per_slice(min_idx:max_idx) = min_band_width;

% repeat with the bottom of k-space
N_bottom = zeros([num_PE/2, num_slices]);
for i = 1:num_slices
    
    for b = 1:num_echos
        % band b
        norm_prob = P_line .* (bands == b);
        norm_prob = norm_prob / sum(norm_prob);
        
        N = pdf_to_mask(norm_prob * shots_per_slice(i), ones(size(norm_prob)));
        
        N_bottom(:,i) = N_bottom(:,i) + N';
    end
   
end


% convert the top mask into shots
PE1 = zeros(sum(shots_per_slice)+sum(shots_per_slice_top), num_echos);
PE2 = zeros(sum(shots_per_slice)+sum(shots_per_slice_top), num_echos);

n = 1;
for i = 1:num_slices

    % do top shots first
    N_i = N_top(:,i);
    r = find(N_i);

    shots = reshape(r, [length(r)/num_echos, num_echos]);
    
    % reverse the order of the shots so that they go from the center
    % outwards, with the shot closest to the center first
    shots = shots(end:-1:1, end:-1:1);
    
    
    numshots = length(r)/num_echos;
    
    PE2(n:n+numshots-1, :) = i;
    PE1(n:n+numshots-1, :) = shots;
    
    PE1(n+numshots-1:-1:n, :) = shots; %IT: Adjust to match fse3d, 23 Nov 14
    
    n = n+numshots;
    
    
    % repeat with bottom shots
    
    N_i = N_bottom(:,i);
    r = find(N_i);

    shots = reshape(r, [length(r)/num_echos, num_echos]);
    
    % reverse the order of the shots so that they go from the center
    % outwards, with the shot closest to the center first
    shots = shots(end:-1:1, end:-1:1);
    
    % for the bottom of k-space, the shots will be going in the other
    % direction
    shots = num_PE/2 + (num_PE/2 - shots) + 1;
    
    numshots = length(r)/num_echos;
    
    PE2(n:n+numshots-1, :) = i;
    PE1(n:n+numshots-1, :) = shots;
    
    n = n+numshots;
    
end


% just as a double check, show each shot in a different colour

mask = zeros(size(P));

for i = 1:size(PE1,1)
    mask(PE1(i,1), PE2(i,1)) = 8;
    mask(PE1(i,2), PE2(i,1)) = 1;
    mask(PE1(i,3), PE2(i,1)) = 5;
    mask(PE1(i,4), PE2(i,1)) = 2;
    mask(PE1(i,5), PE2(i,1)) = 6;
    mask(PE1(i,6), PE2(i,1)) = 3;
    mask(PE1(i,7), PE2(i,1)) = 7;
    mask(PE1(i,8), PE2(i,1)) = 4;
    
end


PE1 = PE1 - num_PE/2;
PE2 = PE2 - num_slices/2;

PE2 = PE2 - 1; %IT: Adjust to match fse3d, 23 Nov 14

