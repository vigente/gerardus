%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Merge multiple surfaces and remove self-intersection elements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% generate a mesh for 
load rat_head.mat
[node0,face0]=v2s(volimage,0.5,2);

c0=mean(meshcentroid(node0,face0(:,1:3)));
c1=2*[18.5 20.7 20.9]-c0;

[cnode,cface]=meshacylinder(c0,c1,4);
[cnode,cface]=meshcheckrepair(cnode,cface(:,1:3));

% combine two surfaces, producing 4 pieces of subsurfaces: surf 1
% outside/inside of surf2 and surf2 outside/inside of surf1

[no,el]=surfboolean(node0,face0(:,1:3),'all',cnode,cface);
figure
plotmesh(no,el,'y>20')

% take the first surface only

% el(:,4)==1: surf 1 outside of surf 2; el(:,4)==3: surf 1 inside of surf 2
[no,el]=surfboolean(node0,face0(:,1:3),'first',cnode,cface);
figure
plotmesh(no,el)

% the mesh after boolean operation can have self-intersecting elements, one
% has to fix those defects before passing to s2m

[no1,el1]=meshcheckrepair(no(:,1:3),el(:,1:3));

%ISO2MESH_TETGENOPT=' -A -q 0.8 -a 10 ';
[node,elem,face]=s2m(no1,el1,1,10);
figure;
plotmesh(node,face)
