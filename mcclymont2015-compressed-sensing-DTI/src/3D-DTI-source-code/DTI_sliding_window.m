function [ K2 ] = DTI_sliding_window(K, bval )
%DTI_SLIDING_WINDOW Fills un-sampled positions in k-space with the nearest
%acquired direction
%   Used by Awate & DiBella, 2013 ISBI
%   Might be good for initialisation

K2 = K;

sz = size(K);

is_3D = false;
if length(sz) == 4
    is_3D = true;
end

for dir = 1:sz(end)
    
    bval_dir = bval(:,:,dir);
    
    bval_diff = bsxfun(@minus, bval, bval_dir);
    
    % sum of squares difference to each of the other directions
    sum_sq = squeeze(sum(sum(bval_diff.^2, 1), 2));
    sum_sq(dir) = inf;
    
    % repeat for the opposite direction (assuming diffusion is symmetric)
    bval_diff = bsxfun(@minus, bval, -bval_dir);
    
    % sum of squares difference to each of the other directions
    sum_sq2 = squeeze(sum(sum(bval_diff.^2, 1), 2));
    sum_sq2(dir) = inf;
    
    sum_sq = min(sum_sq, sum_sq2);
    
    [~, idx] = sort(sum_sq, 'ascend');
    
    for dir2 = idx(1:end-1)' % go to the second last one, because the last is inf
        if is_3D
            % find locations
            locs = (K2(:,:,:,dir) == 0) & (K(:,:,:,dir2) ~= 0);
            % put them in the right place
            Ktemp = K(:,:,:,dir2);
            K2temp = K2(:,:,:,dir);
            K2temp(locs) = Ktemp(locs);
            K2(:,:,:,dir) = K2temp;
        else
            % find locations
            locs = (K2(:,:,dir) == 0) & (K(:,:,dir2) ~= 0);
            % put them in the right place
            Ktemp = K(:,:,dir2);
            K2temp = K2(:,:,dir);
            K2temp(locs) = Ktemp(locs);
            K2(:,:,dir) = K2temp;
        end
    end   
end


end

