%Transforming the images to get exactly the same result we would get after
%performing the same transform in ITK 
function trans_img=apply_ITK_matrices(orig_img, rotITK, transITK, resolution, centered);

if(~centered)
    [x,y,z]=ndgrid(0:resolution(1):(size(orig_img,1)-1)*resolution(1), ...
        0:resolution(2):(size(orig_img,2)-1)*resolution(2), 0:resolution(3):(size(orig_img,3)-1)*resolution(3));
else
    [x,y,z]=ndgrid(-(size(orig_img,1)-1)*resolution(1)/2:resolution(1):(size(orig_img,1)-1)*resolution(1)/2, ...
        -(size(orig_img,2)-1)*resolution(2)/2:resolution(2):(size(orig_img,2)-1)*resolution(2)/2, ...
        -(size(orig_img,3)-1)*resolution(3)/2:resolution(3):(size(orig_img,3)-1)*resolution(3)/2);
end

coords=[reshape(x,[1 prod(size(x))]); reshape(y,[1 prod(size(x))]); reshape(z,[1 prod(size(x))])];
coords_rot=rotITK*coords+repmat(transITK',[1 size(coords,2)]);
xr=reshape(coords_rot(1,:),size(orig_img));
yr=reshape(coords_rot(2,:),size(orig_img));
zr=reshape(coords_rot(3,:),size(orig_img));
clear coords coords_rot
trans_img=interpn(x,y,z,orig_img,xr,yr,zr);
trans_img(isnan(trans_img))=0;

%To obtain the most similar result to ITK, one should also use:
% trans_img = floor(trans_img);
