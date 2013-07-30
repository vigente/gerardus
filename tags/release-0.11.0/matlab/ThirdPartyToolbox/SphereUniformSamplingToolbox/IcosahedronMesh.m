function TR=IcosahedronMesh
% Name speaks for itself.

% Get the vertex coordinates
t=(1+sqrt(5))/2; % golden ratio
x=[0 1 t];
s=[1 1 1; 1 1 -1; 1 -1 -1; 1 -1 1];
x=repmat(x,[4 1]).*s;
x=[x;circshift(x,[0 -1]);circshift(x,[0 -2])];
x_L2=sqrt(sum(x.^2,2));
x=bsxfun(@rdivide,x,x_L2);

% Triangulate the points
Tri = fliplr(convhulln(x));
TR=TriRep(Tri,x);
