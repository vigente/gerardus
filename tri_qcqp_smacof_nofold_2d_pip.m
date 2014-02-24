function [con, bnd] = tri_qcqp_smacof_nofold_2d_pip(tri, ymin, ymax, amin, amax, isFree, y)
% TRI_QCQP_SMACOF_NOFOLD_2D_PIP  Constraints in PIP format for QCQP-SMACOF
% to ensure that a 2D triangular mesh has no fold-overs.
%
% This function generates the linear and quadratic constraints that can be
% passed to qcqp_smacof() to map a triangular mesh onto 2D without
% fold-overs or, equivalently, untangle a 2D projection of a mesh.
%
% [LB, UB, W, A, RL, RU, QC] 
%   = tri_qcqp_smacof_nofold_2d(TRI, YMIN, YMAX, AMIN, AMAX)
%
%   TRI is a 3-column matrix with a surface mesh triangulation. Each row
%   gives the indices of one triangle. The mesh needs to be a 2D manifold,
%   that can be embedded in 2D or 3D space. All triangles must be oriented
%   counter-clockwise (2D space) or with the normals pointing outwards (3D
%   space).
%
%   YMIN, YMAX are 2-vectors that bound the output variables. That is, once
%   the QCQP problem is solved, all the output points will be within a box
%   with bottom-left and top-right coordinates YMIN and YMAX, respectively.
%
%   AMIN, AMAX are scalars with the minimum and maximum area allowed for
%   each output triangle, respectively.
%
%   LB, UB are matrices with the same size as Y. They provide lower and
%   upper bounds, respectively, for the corresponding coordinates.
%
%   W is a weight matrix the same size as D. W_ij = 0 means that the
%   distance between points i and j does not affect the stress measure. All
%   weights between connected nodes are 1.
%
%   A, RL, RU are a matrix and two vectors, respectively, that provide the
%   linear constraints for the SCIP quadratic program objective function, 
%   RL <= A * Y <= RU. In this case, A = [], as this problem does not
%   require linear constraints.
%
%   QC is a struct with the quadratic constraints for the SCIP quadratic
%   program, each of the form qrl <= y' * Q * y + l' * x <= qru. Its fields
%   are:
%
%     'Q': cell vector. Each cell contains the Q matrix for a quadratic
%           constraint.
%     'l': matrix. Each column contains the l vector for a quadratic
%           constraint.
%     'qrl': vector. Each element is the lower bound of a quadratic
%          constraint.
%     'qul': vector. Each element is the upper bound of a quadratic
%          constraint.
%
% ... = tri_qcqp_smacof_nofold_2d(..., ISFREE, Y)
%
%   ISFREE is an boolean N-vector, where N is the number of vertices in the
%   mesh. ISFREE(i)==true means that the i-th vertex is a free vertex (i.e.
%   an unknown in the QCQP problem). ISFREE(i)==false means that the i-th
%   vertex is a fixed vertex (i.e. with known constant coordinates). By
%   default, all vertices are assumed to be free.
%
%   Y is an (N, 2)-matrix that provides the coordinates of the fixed
%   vertices as Y(ISFREE, :). The values Y(~ISFREE, :) are simply ignored.
%   Thus, Y doesn't need to be provided if all vertices are free. If
%   there's at least a fixed vertex, then Y must be provided.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
% Version: 0.0.1
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
% along with this program.  If not, see
% <http://www.gnu.org/licenses/>.

% check arguments
narginchk(5, 7);
nargoutchk(0, 2);

%% Input arguments

if (size(tri, 2) ~= 3)
    error('TRI must have 3 columns')
end
if (~isvector(ymin) || (length(ymin) ~= 2))
    error('YMIN must be a 2-vector')
end
if (~isvector(ymax) || (length(ymax) ~= 2))
    error('YMAX must be a 2-vector')
end
if (~isscalar(amin))
    error('AMIN must be a scalar')
end
if (~isscalar(amax))
    error('AMAX must be a scalar')
end

% number of vertices
N = length(unique(tri(:)));

if (nargin < 7 || isempty(isFree))
    % if the user doesn't specify which vertices are free and which ones
    % are fixed, we assume that all vertices are free
    isFree = true(N, 1);
end

if (~islogical(isFree))
    error('ISFREE must be of type logical (boolean)')
end

% number of triangles
Ntri = size(tri, 1);

% number of free vertices
Nfree = nnz(isFree);

if (any(~isFree))
    
    if (nargin < 7)
        error('If any vertices are non-free, then initial configuration Y must be provided')
    end
    
    % check dimensions of initial configuration
    if ((size(y, 1) ~= N) || (size(y, 2) ~= 2))
        error('User says there is at least one fixed vertex, but initial configuration Y0 either has not been provided or has wrong dimensions')
    end
    
end

% init outputs
con = cell(1, 2*Ntri+1);
con{1} = 'Subject to';
bnd = cell(1, 2*N+1);
bnd{1} = 'Bounds';

%% Weight matrix

% weight matrix (simply, each element is 1 if the two nodes are connected
% in the mesh, and 0 if they are not)
w = dmatrix_mesh(tri);

%% Upper and lower bounds for the objective function variables

% variables bounds, lb<=nu<=ub
for I = 1:N

    % bounds for x-coordinate
    bnd{2*I} = sprintf(...
        ' %.6g <= x%d <= %.6g', ...
        ymin(1), I, ymax(1));

    % bounds for y-coordinate
    bnd{2*I+1} = sprintf(...
        ' %.6g <= y%d <= %.6g', ...
        ymin(2), I, ymax(2));
    
end

%% Quadratic constraints (each triangle produces a constraint)

% Note: we are going to create our matrices and vectors for the whole 2*N
% variables produced by both free and fixed vertices. The reason is that
% this makes indexing a lot simpler. Once the constraints are created, we
% remove the rows and columns that correspond to the fixed vertices (these
% rows and columns have all elements equal to zero)

% each triangle contributes one quadratic constraint of the forms
%   Amin <=            l'*nu + c < Amax     (1 free vertex)
%   Amin <= nu'*Q*nu + l'*nu     < Amax     (2 free vertices)
%   Amin <= nu'*Q*nu             < Amax     (3 free vertices)
for I = 1:Ntri
    
    % coordinates of the three vertices in the triangle
    triloc = tri(I, :);
    isFreeloc = isFree(triloc);
    
    % depending on the number of free vertices in the triangle, we create
    % different constraints
    switch (nnz(isFreeloc))
        
        case 0 % 0 free vertices, 3 fixed vertices
            
            % this case doesn't contribute any constraints to the quadratic
            % program, as the fixed vertices cannot be moved
            qc.Q{I} = sparse(2*N, 2*N); % Q = 0;
            
        case 1 % 1 free vertex, 2 fixed vertices
    
            % position of the free vertex
            idx = find(isFreeloc);
            
            % shift the vertices until the free vertex is the last one,
            % without changing the sign of the area
            yloc = y(triloc, :);
            yloc = circshift(yloc, 3 - idx);
            isFreeloc = circshift(isFreeloc, 3 - idx);
            if any(isFreeloc ~= [0 0 1]')
                error(['Assertion: Triangle ' num2str(I) ...
                    ' vertices cannot be shifted to canonic configuration'])
            end
            triloc = circshift(triloc(:), 3 - idx);
            
            % auxiliary variables to make the code more readable
            xi = yloc(1, 1); % x-coordinate of 1st fixed vertex
            yi = yloc(1, 2); % y-coordinate of 1st fixed vertex
            xj = yloc(2, 1); % x-coordinate of 2nd fixed vertex
            yj = yloc(2, 2); % y-coordinate of 2nd fixed vertex
            
            k = triloc(3);   % index of free vertex
            
            % constraint with lower bound. Example:
            % c1: -2 x1 +3.23 x4 +1 x2 * x3 >= -1
            con{2*I} = sprintf( ...
                ' c%d: %.6g x%d + %.6g y%d + %.6g >= %.6g', ...
                2*I-1, 0.5*(yi-yj), k, 0.5*(xj-xi), k, 0.5*(xi*yj-xj*yi), amin);
            
            % constraint with upper bound. Example:
            % c1: -2 x1 +3.23 x4 +1 x2 * x3 <= 2
            con{2*I+1} = sprintf( ...
                ' c%d: %.6g x%d + %.6g y%d + %.6g <= %.6g', ...
                2*I, 0.5*(yi-yj), k, 0.5*(xj-xi), k, 0.5*(xi*yj-xj*yi), amax);
            
        case 2 % 2 free vertices, 1 fixed vertex
            
            % position of the fixed vertex
            idx = find(~isFreeloc);
            
            % shift the vertices until the fixed vertex is the first one,
            % without changing the sign of the area
            yloc = y(triloc, :);
            yloc = circshift(yloc, 1 - idx);
            isFreeloc = circshift(isFreeloc, 1 - idx);
            if any(isFreeloc ~= [0 1 1]')
                error(['Assertion: Triangle ' num2str(I) ...
                    ' vertices cannot be shifted to canonic configuration'])
            end
            triloc = circshift(triloc(:), 1 - idx);
            
            % auxiliary variables to make the code more readable
            xi = yloc(1, 1); % x-coordinate of fixed vertex
            yi = yloc(1, 2); % y-coordinate of fixed vertex
            
            j = triloc(2);   % index of free vertex
            k = triloc(3);   % index of free vertex
            
            % constraint with lower bound
            con{2*I} = sprintf( ...
                ' c%d: 0.5 x%d y%d - 0.5 x%d y%d - %.6g x%d + %.6g x%d + %.6g y%d - %.6g y%d >= %.6g', ...
                2*I-1, j, k, k, j, 0.5*yi, j, 0.5*yi, k, 0.5*xi, j, 0.5*xi, k, amin);
            
            % constraint with upper bound
            con{2*I+1} = sprintf( ...
                ' c%d: 0.5 x%d y%d - 0.5 x%d y%d - %.6g x%d + %.6g x%d + %.6g y%d - %.6g y%d <= %.6g', ...
                2*I, j, k, k, j, 0.5*yi, j, 0.5*yi, k, 0.5*xi, j, 0.5*xi, k, amax);
            
        case 3 % three free vertices
            
            % auxiliary variables to make the code more readable

            i = triloc(1);   % index of free vertex
            j = triloc(2);   % index of free vertex
            k = triloc(3);   % index of free vertex
            
            % constraint with lower bound
            con{2*I} = sprintf( ...
                ' c%d: -0.5 x%d y%d +0.5 x%d y%d +0.5 x%d y%d -0.5 x%d y%d -0.5 x%d y%d +0.5 x%d y%d >= %.6g', ...
                2*I-1, j, i, k, i, i, j, k, j, i, k, j, k, amin);
            
            % constraint with upper bound
            con{2*I+1} = sprintf( ...
                ' c%d: -0.5 x%d y%d +0.5 x%d y%d +0.5 x%d y%d -0.5 x%d y%d -0.5 x%d y%d +0.5 x%d y%d <= %.6g', ...
                2*I, j, i, k, i, i, j, k, j, i, k, j, k, amax);
            
        otherwise
            
            error('Assertion error: Triangle has more than 3 vertices')
            
    end
    
end

