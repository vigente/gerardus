%MIRT3D_MEXINTERP  Fast 3D linear interpolation 
%  
% Output_image = mirt3D_mexinterp(Input_image, XI,YI,ZI) interpolates the 3D image 'Input_image' at
%    the points with coordinates X,Y,Z. Input_image is assumed to  be defined at a regular grid 1:N, 1:M, 1:K, 
%    where [M,N,K]=size(Input_images). Points outside the boudary return NaNs. 
%    This is equivalent (but much faster) to Matlab's:
%    Output_image = interp3(Input_image,XI,YI,ZI,'linear',NaN);
%  
% Output_images = mirt3D_mexinterp(Input_images, XI,YI,ZI). Input_images can be a stack of many 3D images (4D).
%   The function interpolates each of the 3D images at X,Y,Z coordinates and return a stack of corresponding
%   interpolated images. This is equivalent to Matlab's
% 
%   Input_images(:,:,:,1)=Input_image1;
%   Input_images(:,:,:,2)=Input_image2;
%   Input_images(:,:,:,3)=Input_image3;
%   Input_images(:,:,:,4)=Input_image4;
%  
%   Output_images(:,:,:,1) = interp3(Input_image1,XI,YI,ZI,'linear',NaN);
%   Output_images(:,:,:,2) = interp3(Input_image2,XI,YI,ZI,'linear',NaN);
%   Output_images(:,:,:,3) = interp3(Input_image3,XI,YI,ZI,'linear',NaN);
%   Output_images(:,:,:,4) = interp3(Input_image4,XI,YI,ZI,'linear',NaN);
% 
%  This is especially usefull fpr vector valued 3D images, RGB images, to interpolate the whole 3D video at the same coordinates
%  or to interpolate image and its gradients at the same time (in image registration).
%  The speed gain is also from the precomputation of nearest points for interpolation, which are the same for all images in a stack.
%
%  Andriy Myronenko, Feb 2008, email: myron@csee.ogi.edu, 
%  homepage: http://www.bme.ogi.edu/~myron/


% The function below compiles the mirt3D_mexinterp.cpp file if you haven't done it yet.
% It will be executed only once at the very first run.
function Output_images = mirt3D_mexinterp(Input_images, XI,YI,ZI)

pathtofile=which('mirt3D_mexinterp.cpp');
pathstr = fileparts(pathtofile);
mex(pathtofile,'-outdir',pathstr);

Output_images = mirt3D_mexinterp(Input_images, XI,YI,ZI);

end