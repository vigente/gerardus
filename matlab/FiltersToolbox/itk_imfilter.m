function itk_imfilter
% ITK_IMFILTER: Run ITK filter on a 2D or 3D image
%
% This function is a multiple-purpose wrapper to be able to run all
% ITK filters that inherit from itk::ImageToImageFilter on a Matlab
% 2D image or 3D image volume.
%
% B = ITK_IMFILTER(TYPE, A)
%
%   TYPE is a string with the filter we want to run. Currently, only
%   the following options are implemented:
%
%     'skel': (BinaryThinningImageFilter3D) Skeletonize a
%             segmentation mask
%
%   A is a 2D matrix or 3D volume with the image or
%   segmentation. Currently, A can be of any of the following
%   Matlab classes:
%
%     boolean
%     double
%     single
%     int8
%     uint8
%     int16
%     uint16
%     int32
%     int64
%
%   B has the same size and class as A, and contains the filtered
%   image or segmentation mask.
%
% This function must be compiled before it can be used from Matlab.
% If Gerardus' root directory is e.g. ~/gerardus, type from a
% linux shell
%
%    $ cd ~/gerardus/matlab
%    $ mkdir bin
%    $ cd bin
%    $ cmake ..
%    $ make install
%
% If cmake throws an error because it cannot find Matlab, then edit
% gerardus/matlab/CMakeLists.txt, and where it says
%
%    SET(MATLAB_ROOT "/usr/local/matlab/R2010b/")
%
%  change to your own Matlab root path.
 
% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.1.0
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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

error('MEX file not found')
