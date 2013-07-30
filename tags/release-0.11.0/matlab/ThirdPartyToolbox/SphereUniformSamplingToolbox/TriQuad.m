function TR=TriQuad(TR)
% Subdivide triangular mesh using triangular quadrisection. During this
% operation new vertices are inserted at the edge midpoints thereby 
% producing four new faces for every face of the original mesh. 
% Illustration of this operation is provided below:
% 
%                     x3                        x3
%                    /  \      subdivision     /  \
%                   /    \         -->        v3__v2
%                  /      \                  / \  / \
%                x1________x2              x1___v1___x2
%
%                   Original vertices:    x1, x2, x3
%                   New vertices:         v1, v2, v3
%
% INPUT ARGUMENTS:
%   - TR   : input mesh. TR can be specified as a TriRep object or a cell,
%            such that TR={X Tri}, where X is the list of vertex 
%            coordiantes and Tri is the list of faces.
%
% OUTPUT:
%   - TR  : subdivided mesh. Same format as input.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
% DATE: May.2012
%

% Get the list of vertex cooridnates and list of faces
if strcmp(class(TR),'TriRep')
    X=TR.X;
    Tri=TR.Triangulation;
elseif iscell(TR) && numel(TR)==2
    X=TR{1};
    Tri=TR{2};
    
    % Check the format
    if isempty(X) || ~ismatrix(X) || ~isnumeric(X) || ~strcmp(class(X),'double') || isequal(X,round(X))
        error('Invalid entry for the list of vertices')
    end
    
    if isempty(Tri) || ~ismatrix(X) || ~isnumeric(Tri) || ~isequal(Tri,round(Tri))
        error('Invalid entry for the list of faces')
    end
    
    Tri=double(Tri);
    
else
    error('Unrecognized input format')
end
    
Nx=size(X,1);   % # of vertices
Nt=size(Tri,1); % # of faces

% Compute edge mid-points
V1=(X(Tri(:,1),:)+X(Tri(:,2),:))/2;
V2=(X(Tri(:,2),:)+X(Tri(:,3),:))/2;
V3=(X(Tri(:,3),:)+X(Tri(:,1),:))/2;
V =[V1;V2;V3];

% Remove repeating vertices 
[V,~,idx]=unique(V,'rows');

% Assign indices to the new triangle vertices
V1= Nx + idx(1:Nt);
V2= Nx + idx((Nt+1):2*Nt);
V3= Nx + idx((2*Nt+1):3*Nt);
clear idx

% Define new faces
T1= [Tri(:,1) V1 V3];
T2= [Tri(:,2) V2 V1];
T3= [Tri(:,3) V3 V2];
T4= [V1       V2 V3];
clear V1 V2 V3

% New mesh
X=[X;V]; 
Tri=[T1;T2;T3;T4];

if strcmp(class(TR),'TriRep')
    TR=TriRep(Tri,X);
else
    TR={X,Tri};
end

