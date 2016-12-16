function [cc, vBadIdx] = sphtri_foldcc(tri, x, volTriMin, dcon, isManuallyGood)
% SPHTRI_FOLDCC  Find connected components of groups of vertices causing a
% fold on the sphere
%
% A triangle on the sphere is bad (flipped) if the volume of the
% tetrahedron it forms with the centre of the sphere is < 0. (In practice,
% it's useful to also consider very small positive triangles as bad, see
% VOLTRIMIN below).
%
% A vertex shared by at least one bad triangle is considered a bad vertex.
%
% Additionally, good vertices that are completely surrounded by bad
% vertices are also considered bad vertices. The reason is that a fold that
% folds again can produce good triangles, but in order to deal with the
% fold, all vertices need to be dealt with.
%
% [CC, vBadIdx] = SPHTRI_FOLDCC(TRI, X, VOLTRIMIN,  DCON, ISMANUALLYGOOD)
%
%   TRI is a 3-column matrix. Each row contains the 3 nodes that form one
%   triangular facet in the mesh.
%
%   X is a 3-column matrix. X(i, :) contains the xyz-coordinates of the
%   i-th node in the mesh.
%
%   VOLTRIMIN is the treshold for bad triangles. It can be a scalar or a
%   vector with one element per triangle in TRI. A triangle t is bad if the
%   volume of the tetrahedron formed by the triangle and the centre of the
%   sphere is < VOLTRIMIN. By default, VOLTRIMIN = 0.
%
%   DCON is the adjacency matrix of TRI. If it is not provided, it is
%   computed internally.
%
%   ISMANUALLYGOOD is a vector with vertices that the user wants to label
%   as "good" even if they are shared by any bad flipped triangles. By
%   default, ISMANUALLYGOOD=[] and is ignored.
%
%   CC is a cell vector. Each element contains a list of bad vertices
%   connected to each other.
%
%   vBadIdx are boolean indices of the bad vertices

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2016 University of Oxford
% Version: 0.3.0
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

% check arguments
narginchk(2, 5);
nargoutchk(0, 2);

% defaults
if ((nargin < 3) || isempty(volTriMin))
    volTriMin = 0;
end

% if the user hasn't provided the adjacency matrix, we compute it
if (nargin < 4 || isempty(dcon))
    dcon = dmatrix_mesh(tri);
end
if (nargin < 5)
    isManuallyGood = [];
end

%% partition mesh into good and bad vertices

% number of vertices
N = size(x, 1);

% signed volume of tetrahedra formed by sphere triangles and origin
% of coordinates
volTri = sphtri_signed_vol(tri, x);

% list of triangles that are flipped or too small, bool
% format
triBadIdx = (volTri < volTriMin);

% list of vertices in bad triangles. These are the vertices that need
% untangling, because they are part of the fold (even if they are shared by
% good and bad triangles)
vTriBad = unique(tri(triBadIdx, :));
vBadIdx = false(N, 1);
vBadIdx(vTriBad) = true;

% the user may want to override the automatic detection, and label some
% vertices as good
if (~isempty(isManuallyGood))
    vBadIdx(isManuallyGood) = false;
end

%% remove good vertices within a fold (a fold that folds again may produce
%% good triangles)

% duplicate adjacency matrix and remove the edges connected to bad vertices
dconGood = dcon;
dconGood(vBadIdx, :) = 0;
dconGood(:, vBadIdx) = 0;

% connected components of the good vertices graph
[~, cc] = graphcc(dconGood);


if (~isempty(cc))

    % number of vertices in each component
    len = cellfun(@length, cc);
    
    % we keep only the largest component, assuming that the rest are good
    % vertices within a fold. Change the labelling of those vertices from good
    % to bad
    idx = (len < max(len));
    vWithinFold = cat(1, cc{idx});
    vBadIdx(vWithinFold) = true;
    
end

clear dconGood

%% partition vertices into connected components

% adjacency matrix of bad vertices
dconBad = sparse(size(dcon, 1), size(dcon, 2));
dconBad(vBadIdx, vBadIdx) = dcon(vBadIdx, vBadIdx);

% connected components of the bad vertices graph
[~, cc] = graphcc(dconBad);
