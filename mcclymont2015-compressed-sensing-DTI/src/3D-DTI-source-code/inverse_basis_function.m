function [ Image_space ] = inverse_basis_function( V, AT, sz )

% split V into smaller segments to save ram
% assume the size will always be a multiple of 8

S = zeros(size(V,1), size(AT,2), 'single');

chunksize = size(V,1) / 8;

for i = 1:8
    idx1 = (i-1) * chunksize + 1;
    idx2 = idx1 + chunksize - 1;
    S(idx1:idx2, :) = single( full(V(idx1:idx2, :)) * AT );
end
% equivalent to S = single(full(V) * AT);

Image_space = reshape(S, sz);

end

