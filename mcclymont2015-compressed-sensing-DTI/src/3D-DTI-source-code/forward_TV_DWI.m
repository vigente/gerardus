function [ TV, TV_grad ] = forward_TV_DWI( I )
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

sz = size(I);

if length(sz) == 3 % 2D sequence 
    
    TV_grad = zeros([sz(1), sz(2), 2, sz(3)]);
    
    TV = 0;
    
    for i = 1:sz(3)
        [TV_vol, TV_grad_vol] = forward_TV_2D(I(:,:,i));

        TV = TV + TV_vol;

        TV_grad(:,:,1,i) = single(TV_grad_vol(:,:,1));
        TV_grad(:,:,2,i) = single(TV_grad_vol(:,:,2));
        
    end

    
else % 3D sequence
    

    TV_grad = zeros([sz(1), sz(2), sz(3), 3, sz(4)], 'single');

    TV = 0;

    for i = 1:sz(4)
        [Dx, Dy, Dz, TV_vol] = forward_TV_mex(double(I(:,:,:,i)), 1);

        TV = TV + TV_vol;

        TV_grad(:,:,:,1,i) = single(Dx);
        TV_grad(:,:,:,2,i) = single(Dy);
        TV_grad(:,:,:,3,i) = single(Dz);
    end
end

end