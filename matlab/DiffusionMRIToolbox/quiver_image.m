function [ h ] = quiver_image( IM, VEC, orient, ix )
%UNTITLED3 Summary of this function goes here
%   Input: 3D image, 4D vector field, orientation, index
%   Output: figure handle

if nargin < 3
    orient = 'slice';
end

switch orient
    case{'row','r'}

        dim_to_collapse = 1;

    case{'col', 'c', 'column'}

        dim_to_collapse = 2;

    case{'s', 'slice'}

       dim_to_collapse = 3;

    otherwise         
        error('Orient should be a string, equal to row, r, col, c, column, s, or slice')
end


if nargin < 4

    ix = ceil(size(IM,dim_to_collapse)/2);
    
end
            


% reshape the image and vector field such that the collapseable dimension
% is at the end

dim_order1 = [setdiff(1:3, dim_to_collapse), dim_to_collapse];
IM = permute(IM, dim_order1);
IM = IM(:,:,ix);

dim_order = [setdiff(1:4, dim_to_collapse), dim_to_collapse];
VEC = permute(VEC, dim_order);
VEC = VEC(:,:,:,ix);

[X, Y] = ndgrid(1:size(IM,1), 1:size(IM,2));

im_min = min(IM(:));
im_max = max(IM(:));
im_range = im_max - im_min;


figure
h = imagesc(IM', [im_min - im_range*2, im_max]);
colormap gray
hold on; 
quiver(X,Y, VEC(:,:,dim_order1(1)), VEC(:,:,dim_order1(2)), 'k');

% to symmetrise:
% quiver(X,Y, -VEC(:,:,dim_order1(1)), -VEC(:,:,dim_order1(2)), 'k');


end

