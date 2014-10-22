function [im_out,rows,cols]=imagesc_nD(im_in,cols,clim)

% Reformats nD-arrayed 2D image data for display as 2D montage.
%
% [IM_OUT,ROWS,COLS]=imagesc_nD(IM_IN,COLS,CLIM)
% 
% IM_IN  = nD data where first two dimensions correspond to a 2D image. 3rd to nth dimension are arrays of 2D images.
% COLS   = Specify number of columns in output image (Optional)
% CLIM   = [CLOW CHIGH] to specify scaling (Optional)
%  
% IM_OUT = 2D mosaic image
% ROWS   = Number of rows in output image
% COLS   = Number of columns in output image
%
% USAGE
% load mri
% imagesc_nD(D);

% Author: Irvin Teh <irvin@well.ox.ac.uk>
% Copyright © 2010-2014 University of Oxford
% Version: 0.1.0
% $Rev$
% $Date$
% 
% University of Oxford means the Chancellor, Masters and Scholars of
% the University of Oxford, having an administrative office at
% Wellington Square, Oxford OX1 2JD, UK. 
%
% This file is part of Gerardus.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. The offer of this
% program under the terms of the License is subject to the License
% being interpreted in accordance with English Law and subject to any
% action against the University of Oxford being under the jurisdiction
% of the English Courts.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% _________________________________________________________________

% Check arguments
narginchk(1, 3);
nargoutchk(0, 3);

% Reshape nD data into 3D data
siz=size(im_in);
x=siz(1);
y=siz(2);
num_images=prod(siz(3:end));
im_in=reshape(im_in,[x y num_images]);

% Determine rows and columns based on screen size
if nargin<2
    siz_screen=get(0,'ScreenSize');
    aspect_ratio=siz_screen(3)/siz_screen(4); %width/height
    cols=ceil(sqrt(x/y*num_images*aspect_ratio));
end
rows=ceil(num_images/cols);

% Array multislice data
im_out=zeros(x*rows,y*cols);
count=0;
for i=1:rows
    for j=1:cols
        count=count+1;
        if count<=num_images
            im_out((i-1)*x+(1:x),(j-1)*y+(1:y))=im_in(:,:,count);
        end
    end
end

% Display images
if nargin<3
    imagesc(im_out)
else
    imagesc(im_out,clim);
end
axis image

end

