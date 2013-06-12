function [xsax, d] = thickslice_collate_sax_la(xsax, d, xla, K)
% THICKSLICE_COLLATE_SAX_LA  Build a distance/adjacency matrix collating
% points from a Long Axis plane to a set of Short Axis planes
%
% [XOUT, DOUT] = thickslice_collate_sax_la(XSAX, DSAX, XLA, K)
%
%   XSAX is a 3-column matrix where each row has the coordinates of a set
%   of points distributed on short axis (SAX) slices.
%
%   DSAX is a square matrix with the distance/adjacency between the points
%   in XSAX.
%
%   XLA is a 3-column matrix with the coordinates of a set of points
%   distributed on a long axis (LA) slice. The LA slice is usually assumed
%   to be more or less orthogonal to the SAX slices, although it doesn't
%   need to be exactly so.
%
%   K is a scalar with the number of nearest neighbours that are
%   considered. By default, K=3.
%
%   XOUT is a 3-colum matrix with the concatenation of XSAX and XLA.
%
%   DOUT is the distance/adjacency matrix for the collated SAX and LA
%   points.
%
% See also: scimat_dmatrix_thickslice

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2013 University of Oxford
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
narginchk(3, 4);
nargoutchk(0, 2);

if (size(xsax, 2) ~= 3 || size(xla, 2) ~= 3)
    error('XSAX and XLA must be 3-column matrices')
end

% number of SAX points
Nsax = size(xsax, 1);

if (Nsax ~= size(d, 1) || Nsax ~= size(d, 1))
    error('DSAX must be a square matrix with length equal to the number of rows in XSAX')
end

% defaults
if (nargin < 4 || isempty(K))
    K = 3;
end

% number of LA points
Nla = size(xla, 1);

% enlarge distance matrix to accommodate the LA points
d = [d zeros(Nsax, Nla); zeros(Nla, Nsax) zeros(Nla, Nla)];

% compute distance matrix within the LA points
ds = dmatrix(xla');

% compute local neighbourhood for the slice and update the big distance
% matrix with it
d = update_dmatrix(d, ds, (Nsax+1:Nsax+Nla)', (Nsax+1:Nsax+Nla)', K);

% compute distance matrix from LA to SAX points
ds = dmatrix(xla', xsax');

% compute local neighbourhood between LA and SAX and update the big
% distance matrix with it
d = update_dmatrix(d, ds, (Nsax+1:Nsax+Nla)', (1:Nsax)', K);

% collate the points themselves
xsax = [xsax; xla];

end

% this function also in scimat_dmatrix_thickslice()
function d = update_dmatrix(d, ds, idx1, idx2, K)

if (size(ds, 1) ~= length(idx1) || size(ds, 2) ~= length(idx2))
    error('Internal error: Slice distance matrix has incongruent dimensions with the index vectors')
end

% sort distances from smaller to larger in each row
[ds, idxto] = sort(ds, 2, 'ascend');

% the first column will be distance 0 of each point to itself
ds(:, 1) = [];
idxto(:, 1) = [];

% transfer only the distances and connections of the K-nearest
% neighbours to the global distance matrix (making sure that the output
% distance matrix is symmetric)
idxto = idxto(:, 1:K);
idxto = idx2(idxto);
idxfrom = repmat(idx1, 1, K);
d(sub2ind(size(d), idxfrom, idxto)) = ds(:, 1:K);
d(sub2ind(size(d), idxto, idxfrom)) = ds(:, 1:K);

end
