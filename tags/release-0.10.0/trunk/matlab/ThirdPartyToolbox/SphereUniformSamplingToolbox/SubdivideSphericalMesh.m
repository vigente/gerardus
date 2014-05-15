function TR=SubdivideSphericalMesh(TR,k)
% Subdivide triangular mesh representing the surface of the unit sphere
% k times using triangular quadrisection (see function TriQuad for more 
% info). The newly inserted vertices are re-projected onto the unit sphere
% after every iteration.
%
% INPUT ARGUMENTS:
%   - TR   : input mesh. TR can must be specified as a TriRep object.
%   - k    : desired number of subdivisions. k=0 is default.
%
% OUTPUT:
%   - TR  : subdivided mesh. Same format as input.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% DATE: June.2012
%

if nargin<2 || isempty(k), k=0; end

if k==0, return; end
for i=1:k
    
    % Subdivide the mesh
    TR=TriQuad(TR);
    
    % Project the points onto the surface of the unit sphere
    x=TR.X;
    x_L2=sqrt(sum(x.^2,2));
    x=bsxfun(@rdivide,x,x_L2);
    TR=TriRep(TR.Triangulation,x);
    
end

