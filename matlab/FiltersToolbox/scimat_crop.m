function scimat = scimat_crop(scimat, from, to, ds)
% SCIMAT_CROP  Crop a SCIMAT image or segmentation volume of any dimension.
%
% This function has two modes: Either crop a 1D to 5D scimat volume, or
% crop a 3D volume given slice-by-slice by a list of image files. The
% latter is convenient for creating virtual slices from a stack of high
% resolution histology images.
%
% -------------------------------------------------------------------------
% Syntax to crop a 1D to 5D scimat volume:
%
% SCIMAT2 = scimat_crop(SCIMAT, FROM, TO)
%
%   SCIMAT is a struct with a metainfo-enriched image or segmentation (see
%   "help scimat" for details). The image can have any number of
%   dimensions.
%
%   FROM, TO are (row, col, slice) vectors with the index coordinates that
%   define the cropping box. The number of indices depend on the
%   dimensionality of the image. For example, for a 3D volume, FROM=[2 3
%   7], TO=[15 20 22].
%
%   If FROM or TO don't have all the dimensions, those dimensions are not
%   cropped. E.g. in an image with size [256 512 128 10],
%
%     FROM=[50 75];
%     TO=[100 89];
%
%   is the same as
%
%     FROM=[50 75 1 1];
%     TO=[100 89 128 10];
%
%   If FROM=[] and TO=[], or they are not given as input argument, then the
%   image is not cropped.
%
%   SCIMAT2 is the output cropped volume.
%
% Example:
%
%   % create example scimat struct
%   scimat = scimat_im2scimat(zeros(256, 512, 128), [2 3 5], [-3.0 7.0 11.0]);
%
%   % crop struct
%   scimat2 = scimat_crop(scimat, [50 75], [100 89]);
%
% -------------------------------------------------------------------------
% Syntax to crop a 3D volume given slice-by-slice by a list of image files:
%
% SCIMAT = scimat_crop(FILE, FROM, TO, DS)
%
%   FILE is a cell vector with a list of path/filenames. Each filename
%   corresponds to a 2D slice of the 3D volume. The files can be any format
%   that scimat_load() can read.
%
%   FROM/TO are 3-vectors with the (row, col, slice) indices of the
%   cropping. NaN values can be provided to mean first/last voxel.
%
%   FROM/TO can also be given as matrices. In that case, each row
%   corresponds to a different crop of the 3D volume. That is, SCIMAT(i)
%   corresponds to FROM(i,:), TO(i,:).
%
%   DS is a scalar with the slice thickness. By default, DS=1.

% Authors: Ramon Casero <rcasero@gmail.com>, Benjamin Villard <b.016434@gmail.com>
% Copyright © 2011-2016 University of Oxford
% Version: 0.5.2
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

switch (class(scimat))
    
    case 'struct' % scimat structure
        
        % check arguments
        narginchk(1, 3);
        nargoutchk(0, 1);
        
        % defaults
        if (nargin < 2)
            from = [];
        end
        if (nargin < 3)
            to = [];
        end
        
        % FROM and TO can be shorter than the number of dimensions of the image,
        % but not longer
        if (length(from) > ndims(scimat.data))
            error('FROM has more dimensions than the image')
        end
        if (length(to) > ndims(scimat.data))
            error('TO has more dimensions than the image')
        end
        
        % default: extend "from" and "to" with dummy dimensions, if necessary, so
        % that they are the same dimension as the image
        % E.g. in an image with size [256 512 128 10],
        %
        %     FROM=[50 75];
        %     TO=[100 89];
        %
        %   is the same as
        %
        %     FROM=[50 75 1 1];
        %     TO=[100 89 128 10];
        sz = size(scimat.data);
        sz(length(sz)+1:ndims(scimat.data)) = 1;
        from(length(from)+1:ndims(scimat.data)) = 1;
        to(length(to)+1:ndims(scimat.data)) = sz(length(to)+1:end);
        
        % "bottom-left" coordinates of what is going to become the first voxel of
        % the volume after cropping
        %
        % swap X and Y coordinates, so that they match the scimat convention for
        % the axes: (r, c, s, f) <-> (y, x, y, t)
        xmin = scimat_index2world(from, scimat);
        xmin = xmin([2 1 3:end]) - [scimat.axis.spacing]/2;
        
        % crop and correct metainformation in the scimat image
        for I = 1:length(scimat.axis)
            
            % crop one dimension of the image. We are going to be shifting
            % circularly the dimensions of the image so that we always crop the
            % first dimension
            idx = cell(1, length(scimat.axis));
            idx{1} = from(I):to(I);
            idx(2:end) = {':'};
            scimat.data = scimat.data(idx{:});
            
            % size
            scimat.axis(I).size = size(scimat.data, 1);
            
            % "left" edge of first voxel
            scimat.axis(I).min = xmin(I);
            
            % circular shift of the dimensions of image so that we can crop along
            % the first dimension in the next iteration
            idx = 1:length(scimat.axis);
            scimat.data = permute(scimat.data, idx([2:end 1]));
            
        end
        
        % Note: the last circular shift of the image dimension has already returned
        % the image to its original shape, so no need to do anything after the loop

    case 'cell' % list of files
        
        % check arguments
        narginchk(1, 4);
        nargoutchk(0, 1);
        
        if (nargin < 4 || isempty(ds))
            ds = 1;
        end
        
        % change nomenclature
        file = scimat;
        clear scimat
        
        % number of files
        NS = length(file);
        
        if (NS == 0)
            scimat = [];
            return
        end
        
        % load first file, to know the image size
        scimat_in = scimat_load(file{1});
        NR = size(scimat_in.data, 1);
        NC = size(scimat_in.data, 2);
        NCh = size(scimat_in.data, 5);
        
        if (size(from, 1) ~= size(to, 1))
            error('FROM and TO must have the same number of rows, one from per crop')
        end
        if (size(from, 2) ~= 3 || size(to, 2) ~= 3)
            error('FROM and TO must have 3 columns')
        end
        
        % replace slice NaN indices by whole range of slices
        from(isnan(from)) = 1;
        to(isnan(to(:, 1)), 1) = NR;
        to(isnan(to(:, 2)), 2) = NC;
        to(isnan(to(:, 3)), 3) = NS;
        
        % box that contains all the crops
        fromMin = min(from);
        toMax = max(to);
        
        % prepare each of the crops (e.g. filling in missing parameters)
        offset = zeros(size(from, 1), 3);
        for N = 1:size(from, 1)
            
            % offset of cropped volume
            offset(N, 1:2) = scimat_index2world(from(N, 1:2), scimat_in);
            offset(N, 3) = (from(N, 3) - 1) * ds;
            
            % init output
            scimat(N) = scimat_im2scimat(...
                zeros(...
                to(N, 1)-from(N, 1)+1, ...
                to(N, 2)-from(N, 2)+1, ...
                to(N, 3)-from(N, 3)+1, 1, NCh, 'like', scimat_in.data), ...
                [scimat_in.axis.spacing ds], offset(N, :));
            
        end
        
        % loop files
        for I = fromMin(3):toMax(3)

            % load file with the 2D slice
            scimat_in = scimat_load(file{I});
            
            % loop crops
            for N = 1:size(from, 1)
                
                % transfer cropped slice
                scimat(N).data(:, :, I-from(N, 3)+1, :, :) ...
                    = scimat_in.data(from(N, 1):to(N, 1), ...
                    from(N, 2):to(N, 2), :, :, :);
                
            end
            
        end

    otherwise
        
        error('First input must be a SCIMAT volume or a CELL vector with a list of filenames')
        
end
