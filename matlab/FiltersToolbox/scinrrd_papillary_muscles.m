function nrrd2 = scinrrd_papillary_muscles(nrrd, NPAPS)
% SCINRRD_PAPILLARY_MUSCLES  Extract the papillay muscles from a
% segmentation of the Left Ventricle's cavity
%
% NRRD = SCINRRD_PAPILLARY_MUSCLES(NRRD)
%
%   This function extracts a segmentation of the papillary muscles from the
%   segmentation of a Left Ventricular cavity. In fact, the papillary
%   muscles are extraced, plus some bits of the wall that don't belong to
%   them, and instead of stopping at the chordae tendineae, the
%   segmentation continues and selects some trabeculations.
%
%   NRRD is the SCI NRRD struct with the segmentation.
%
% NRRD = SCINRRD_PAPILLARY_MUSCLES(NRRD, NPAPS)
%
%   NPAPS is a constant with the number of papillary muscles we want to
%   extract. By default, NPAPS=2.
%
%   Note on SCI NRRD: Software applications developed at the University of
%   Utah Scientific Computing and Imaging (SCI) Institute, e.g. Seg3D,
%   internally use NRRD volumes to store medical data.
%
%   When label volumes (segmentation masks) are saved to a Matlab file
%   (.mat), they use a struct called "scirunnrrd" to store all the NRRD
%   information:
%
%   >>  scirunnrrd
%
%   scirunnrrd = 
%
%          data: [4-D uint8]
%          axis: [4x1 struct]
%      property: []
%

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2010 University of Oxford
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

% check arguments
error( nargchk( 1, 2, nargin, 'struct' ) );
error( nargoutchk( 0, 1, nargout, 'struct' ) );

% defaults
if (nargin < 2 || isempty(NPAPS))
    NPAPS = 2;
end

% remove the dummy dimension
nrrd = scinrrd_squeeze( nrrd );

% get a tight box around the LV segmentation
box = scinrrd_box(nrrd);

% convert real world coordinates to indices
boxi = scinrrd_world2index( box', nrrd.axis )';

% get data volume size
sz = size( nrrd.data );

% create output segmentation volume
nrrd2 = nrrd;
nrrd2.data(:) = 0;

%% Compute convex hull for the segmentation in each slice, and remove the
%% segmentation. This way, we get the papillary muscles together with some
%% rubbish where the convex hull fills in incorrectly a region outside the
%% cavity edge

% element for morphological operator
se1 = strel('disk', 1);

% loop every slice with at least one LV cavity voxel
for I = min(boxi(3,:)):max(boxi(3,:))
    
    % connected components filter (nlabs is the number of labels w/o counting
    % the background)
    [ im, nlabs ] = bwlabel( nrrd.data(:, :, I) );
    
    % get number of voxels in each label
    nvox = zeros( nlabs, 1 );
    for L = 1:nlabs
        nvox(L) = sum( im(:) == L );
    end
    
    % we assume that the largest connected component is the cavity
    [ nvoxmax, idx ] = max( nvox );
    
    % it only makes sense to compute the convex hull if there are at least
    % three pixels
    if ( nvoxmax >= 3 )
    
        % we delete the rest of voxels; even if they are part of the cavity, we
        % don't want them interfering with the convex hull
        im = (im == idx);
        
        % get linear indices of cavity pixels
        idx = find( im );
        
        % convert linear index to multiple subscripts
        [ir, ic] = ind2sub( sz(1:2), idx );
        
        % compute the convex hull
        idxhull = convhull( ir, ic );
        
        % get all the pixels inside of convex hull
        hull = roipoly(im, ic(idxhull), ir(idxhull));
        
        % we are going to erode the hull polygon 1 pixel, to help in some
        % cases to avoid that the fringe of the next XOR operation are
        % connected between them
        hull = imerode(hull, se1);
        
        % put the new segmentation into the data volume
        nrrd2.data(:, :, I) = bitxor( im, hull );
        
    end
    
end

%% Now we are going to get a guess at which bits are papillary muscle and
%% which bits are rubbish. For this, we start start at the middle slice,
%% because this one is likely to show the papillary muscles well. We assume
%% that the largest connected components belong to the papillary muscles,
%% and then look for the components that are touching the previous ones in
%% the next slice

% extract the middle slice
im = nrrd2.data(:, :, round(mean(boxi(3,:))));

% compute connected components
[im, nlabs] = bwlabel( im );

% get number of voxels in each label
nvox = zeros( nlabs, 1 );
for L = 1:nlabs
    nvox(L) = sum( im(:) == L );
end
    
% order connected components by decreasing area size
[nvox, idxlabs] = sort(nvox, 1, 'descend');

% the number of papillary muscles at the middle of the cavity may reduce as
% we move towards the end of the volume, and the muscles fuse together
npaps = NPAPS;

% keep only a number of largest components, each one assumed to belong to a
% papillary muscle
im2 = im*0;
for I = 1:npaps
    im2(im == idxlabs(I)) = I;
end

for I = round(mean(boxi(3,:)))+1:max(boxi(3,:))
   
    % previous slice
    im1 = im2;
    
    % initialize mask that we are going to use to discard pixels that are
    % far from the pixels in the previous slice
    im1mask = im1*0;
    
    % dilate each previous papillary muscle component to create a mask for
    % the region of interest
    for L = 1:npaps
        
        % isolate current component
        im1comp = im1 == L;
        
        % number of voxels in the component
        nvox = length(find(im1comp));
        
        % assuming a quasi circular shape, then the area is a=pi*r^2. Because
        % the area is equal to the number of pixels, we can compute the radius
        % of the component as r=sqrt(nvox/pi)
        rad = sqrt(nvox / pi);
        
        % we are going to allow the component to search for a new component
        % within an area that is the 25% dilation of the old component; but
        % first we do an erosion of 25% to remove spurious bits of
        % segmentation
        im1comp = imerode(im1comp, strel('disk', round(rad*0.25)));

        % dilate the previous component
        im1comp = imdilate(im1comp, strel('disk', round(rad*0.50)));
        
        % add dilated component to the mask
        im1mask = im1mask | im1comp;
    end
    
    % remove all pixels from the new slice that are outside of the mask
    nrrd2.data(:, :, I) = nrrd2.data(:, :, I) & im1mask;
    
    % compute connected components in the new slice
    [nrrd2.data(:, :, I), nlabs2] = bwlabel( nrrd2.data(:, :, I) );
    
    % we need a temporal image to store the new papillary bits
    im2 = nrrd2.data(:, :, I) * 0;

    % make sure that we are not in a slice without any segmented pixels
    if (nlabs2)
        
        % matrix to keep track of the intersection of each papillary bit
        % with every connected component
        sim = zeros( npaps, nlabs2 );
        
        % intersect each papillary bit with each connected component bit
        for L = 1:npaps
            % compute how many pixels of the new component overlap the old
            % component
            for C = 1:nlabs2
                sim(L, C) = length(find((im1 == L) ...
                    & (nrrd2.data(:, :, I) == C)));
            end
        end
        
        % for each papillary, select the connected component with the
        % largest intersection with the previous slice's bit
        for L = 1:npaps
            [foo, Lnew] = max(sim(L, :));
            im2(nrrd2.data(:, :, I) == Lnew) = L;
        end
        
    end
    
    % overrite the segmentation volume with just the papillary bits
    nrrd2.data(:, :, I) = im2;
    
    % recheck the number of components now in the slice (because the
    % papillary muscles can fuse together)
    [ foo, npaps ] = bwlabel( nrrd2.data(:, :, I) );
        
end


% keep only a number of largest components, each one assumed to belong to a
% papillary muscle
im2 = im*0;
for I = 1:NPAPS
    im2(im == idxlabs(I)) = I;
end

% the number of papillary muscles at the middle of the cavity may reduce as
% we move towards the end of the volume, and the muscles fuse together
npaps = NPAPS;

for I = round(mean(boxi(3,:))):-1:min(boxi(3,:))
   
    % previous slice
    im1 = im2;
    
    % initialize mask that we are going to use to discard pixels that are
    % far from the pixels in the previous slice
    im1mask = im1*0;
    
    % dilate each previous papillary muscle component to create a mask for
    % the region of interest
    for L = 1:npaps
        
        % isolate current component
        im1comp = im1 == L;
        
        % number of voxels in the component
        nvox = length(find(im1comp));
        
        % assuming a quasi circular shape, then the area is a=pi*r^2. Because
        % the area is equal to the number of pixels, we can compute the radius
        % of the component as r=sqrt(nvox/pi)
        rad = sqrt(nvox / pi);
        
        % we are going to allow the component to search for a new component
        % within an area that is the 25% dilation of the old component; but
        % first we do an erosion of 25% to remove spurious bits of
        % segmentation
        im1comp = imerode(im1comp, strel('disk', round(rad*0.25)));

        % dilate the previous component
        im1comp = imdilate(im1comp, strel('disk', round(rad*0.50)));
        
        % add dilated component to the mask
        im1mask = im1mask | im1comp;
    end
    
    % remove all pixels from the new slice that are outside of the mask
    nrrd2.data(:, :, I) = nrrd2.data(:, :, I) & im1mask;
    
    % compute connected components in the new slice
    [nrrd2.data(:, :, I), nlabs2] = bwlabel( nrrd2.data(:, :, I) );
    
    % we need a temporal image to store the new papillary bits
    im2 = nrrd2.data(:, :, I) * 0;

    % make sure that we are not in a slice without any segmented pixels
    if (nlabs2)
        
        % matrix to keep track of the intersection of each papillary bit
        % with every connected component
        sim = zeros( npaps, nlabs2 );
        
        % intersect each papillary bit with each connected component bit
        for L = 1:npaps
            % compute how many pixels of the new component overlap the old
            % component
            for C = 1:nlabs2
                sim(L, C) = length(find((im1 == L) ...
                    & (nrrd2.data(:, :, I) == C)));
            end
        end
        
        % for each papillary, select the connected component with the
        % largest intersection with the previous slice's bit
        for L = 1:npaps
            [foo, Lnew] = max(sim(L, :));
            im2(nrrd2.data(:, :, I) == Lnew) = L;
        end
        
    end
    
    % overrite the segmentation volume with just the papillary bits
    nrrd2.data(:, :, I) = im2;
    
    % recheck the number of components now in the slice (because the
    % papillary muscles can fuse together)
    [ foo, npaps ] = bwlabel( nrrd2.data(:, :, I) );
        
end

