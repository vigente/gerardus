function nrrdsk = scinrrd_skeleton_prune(nrrdsk, nrrd, maxclump, minlen, lratio, alphamax, p)
% SCINRRD_SKELETON_PRUNE  Prune branches in a segmentation's skeletonization
%
% This function prunes the leaves of a segmentation skeleton in three
% steps:
%
%   1. Remove clumps of voxels
%
%   2. Prune very short leaves (iteratively until no more leaves can be
%      pruned)
%
%   3. Prune leaves created by artifacts in the segmentation. Spurious
%      leaves are those roughly as long as the local radius of the main
%      branch they are attached to (iteratively until no more leaves can
%      be pruned)
%
%
% NRRDPR = SCINRRD_SKELETON_PRUNE(NRRDSK)
% NRRDPR = SCINRRD_SKELETON_PRUNE(NRRDSK, [], MAXCLUMP, MINLEN)
%
%   This syntax runs steps 1 and 2 only.
%
%   NRRDSK is an SCI NRRD struct. NRRDSK.data contains the result of
%   running a skeletonization algorithm on a binary segmentation NRRD,
%   e.g.
%
%     >> nrrdsk = nrrd;
%     >> nrrdsk.data = itk_imfilter('skel', nrrd);
%
%   MINLEN is a scalar with the minimum length in voxels for a leaf. Any
%   leaf shorter than MINLEN will be pruned. By default, MINLEN = 5 voxel.
%
%
% NRRDPR = SCINRRD_SKELETON_PRUNE(NRRDSK, NRRD, MAXCLUMP, MINLEN, LRATIO, ALPHAMAX, P)
%
%   This syntax runs steps 1, 2 and 3.
%
%   NRRD is the binary segmentation mentioned above.
%
%   LRATIO is a scalar. Leaves with L/R < LRATIO will be pruned, where
%   L is the length from the bifurcation to the tip of the leaf. By
%   default, LRATIO=1.2.
%
%   ALPHAMAX and P are merging parameters. See the help of function
%   skeleton_label for details. If merging is not enabled, then no leaves
%   will be pruned in step 3.
%
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

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2011 University of Oxford
% Version: 0.3.0
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
error(nargchk(1, 7, nargin, 'struct'));
error(nargoutchk(0, 1, nargout, 'struct'));

if ~isstruct(nrrdsk)
    error('NRRDSK must be an SCI NRRD struct')
end

% defaults
if (nargin < 2)
    nrrd = [];
end
if (nargin < 3 || isempty(maxclump))
    maxclump = 9;
end
if (nargin < 4 || isempty(minlen))
    minlen = 5;
end
if (nargin < 5 || isempty(lratio))
    lratio = 1.2;
end
if (nargin < 6 || isempty(alphamax))
    alphamax = -Inf; % don't merge by default
end
if (nargin < 7 || isempty(p))
    p = 1.0; % don't smooth by default
end

%% Step 1: removal of big clumps of voxels

% label bifurcation clumps
[~, ~, bifcc] = skeleton_label(nrrdsk.data, [], [nrrdsk.axis.spacing]);

% get number of voxels in each clump
n = cellfun(@(x) length(x), bifcc.PixelIdxList);

% find clumps larger than maxclump
idx = find(n > maxclump);

% remove clumps larger than maxclump
nrrdsk.data(cat(1, bifcc.PixelIdxList{idx})) = 0;

%% Step 2: pruning of very short leaf branches
while (1)
    
    % compute skeleton labelling
    [~, cc] = skeleton_label(nrrdsk.data, [], [nrrdsk.axis.spacing]);
   
    % get number of voxels in each branch
    n = cellfun(@(x) length(x), cc.PixelIdxList);
    
    % find leaf-branches that are shorter than the minimum length
    idx1 = find(n < minlen & cc.IsLeaf);
    
    % remove short branches from the segmentation
    nrrdsk.data(cat(1, cc.PixelIdxList{idx1})) = 0;
    
    % recompute the skeleton labelling
    [~, ~, bifcc, mcon] = skeleton_label(nrrdsk.data, [], [nrrdsk.axis.spacing]);
   
    % find bifurcation clusters that are connected to 0 or 1 branches
    idx2 = find(sum(mcon, 1) < 2);
    
    % remove those bifurcation clumps, because they are not connecting
    % branches, they are either floating alone in space, or terminating a
    % branch
    nrrdsk.data(cat(1, bifcc.PixelIdxList{idx2})) = 0;
    
    % if no bifurcation clumps were found to be removed, stop the
    % algorithm, because that means that no new short leaf-braches can be
    % found either
    if (isempty(idx2))
        break;
    end
    
end


%% Step 3: pruning of leaf branches produced by segmentation artifacts

% skip if no full segmentation is provided
if isempty(nrrd)
    return
end

% repeat the process until no branches are removed
atleastonepruning = true;
while (atleastonepruning)
    disp('hi')
    
    % if there are no prunings in this iterations, we stop
    atleastonepruning = false;
    
    % label the segmentation using multiple merging at every bifurcation
    % clump
    [nrrd.data, cc, bifcc, mcon, madj, cc2, mmerge] = ...
        skeleton_label(nrrdsk.data, nrrd.data, [nrrd.axis.spacing], ...
        alphamax, p, false);
    
    % measure the stats of every merged branch
    stats2 = scinrrd_seg2label_stats(nrrd, cc2, p);
    
    % loop each branch in the list of merged branches
    for I = 1:cc2.NumObjects
        
        % we consider the current branch a main branch, and any branch
        % coming out of it, a secondary branch
        
        % compute the major radius of the main branch
        r = sqrt(4 * stats2.var(2, I));
        
        % loop consecutive pairs of segments in the main branch
        for J = 1:length(cc2.MergedBranches{I})-1
            
            % get the bifurcation clump between both segments
            bif = cc2.MergedBifClumps{I}(J);
            
            % get secondary branches attached to the bifurcation clump
            % (usually, there's only 1 secondary branch, but there can be more)
            bn = setdiff(find(mcon(:, bif)), cc2.MergedBranches{I}(J:J+1));
            
            % keep only secondary branches that are leaf-branches
            bn = bn(cc.IsLeaf(bn));
            
            % skip to next segment if there are no valid leaf secondary
            % branches
            if isempty(bn)
                continue
            end
            
            % compute distance between the first and last voxels in the
            % secondary branch, which gives a better estimate of whether the
            % branch protudes from the main vessel or not, than the branches
            % length (because the branch can be bent)
            len = zeros(size(bn));
            for K = 1:length(len)
                idx = cc.PixelIdxList{bn(K)}([1 end]);
                [r1, c1, s1] = ind2sub(size(nrrdsk.data), idx);
                xyz = scinrrd_index2world([r1, c1, s1], nrrdsk.axis);
                len(K) = sqrt(sum(diff(xyz).^2));
            end
            
            % get branches that have to be pruned
            toremove = (len / r) <= lratio;
            atleastonepruning = atleastonepruning || any(toremove);
            
            % if length of the secondary branch is not much larger than the
            % radius of the main branch, we assume that the secondary branch is
            % a segmentation artifact, and remove it
            nrrdsk.data(cat(1, ...
                cc.PixelIdxList{bn(toremove)})) = 0;
            
        end
        
    end
    
end
