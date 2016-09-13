% qmeshcut demonstration
%
% by Qianqian Fang, <fangq at nmr.mgh.harvard.edu>
%
% to demonstrate how to use qmeshcut to produce cross-sectional plot 
% of an un-structured (tetrahedral) mesh

% run vol2mesh demo 1 to create a 3d mesh

demo_vol2mesh_ex1

% define a plane by 3 points, in this case, z=mean(node(:,3))

z0=mean(node(:,3));

plane=[min(node(:,1)) min(node(:,2)) z0
       min(node(:,1)) max(node(:,2)) z0
       max(node(:,1)) min(node(:,2)) z0];

% run qmeshcut to get the cross-section information at z=mean(node(:,1))
% use the x-coordinates as the nodal values

[cutpos,cutvalue,facedata]=qmeshcut(elem(:,1:4),node,node(:,1),plane);

% plot your results

figure;
hsurf=trimesh(face(:,1:3),node(:,1),node(:,2),node(:,3),'facecolor','none');
hold on;
if(isoctavemesh)
  hcut=patch('Faces',facedata,'Vertices',cutpos);
else
  hcut=patch('Faces',facedata,'Vertices',cutpos,'FaceVertexCData',cutvalue,'facecolor','interp');
end
%set(hcut, 'linestyle','none')
axis equal;

% qmeshcut can also cut a surface

[bcutpos,bcutvalue,bcutedges]=qmeshcut(face(:,1:3),node,node(:,1),plane);
[bcutpos,bcutedges]=removedupnodes(bcutpos,bcutedges);
bcutloop=extractloops(bcutedges);

bcutloop(isnan(bcutloop))=[]; % there can be multiple loops, remove the separators

% plot the plane-surface cuts

plot3(bcutpos(bcutloop,1),bcutpos(bcutloop,2),bcutpos(bcutloop,3),'r','LineWidth',4);

% essencially, this should be the same as you do a removedupnodes(cutpos,facedata)
% and then call extractloop(facedata)


% qmeshcut can also cut along an isosurface

% define a field over the mesh: sensitivity map from a source/detector pair

r1=[node(:,1)-20,node(:,2)-25,node(:,3)-25];
r2=[node(:,1)-10,node(:,2)-25,node(:,3)-14];
r1=sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2);
r2=sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2);

k=10;
g1=exp(sqrt(-1)*k*r1)./(4*pi*r1); % calculate the Green's function
g2=exp(sqrt(-1)*k*r2)./(4*pi*r2);
g12=g1.*g2;  % this is the sensitivity map

figure
plotmesh([node log10(abs(g12))],elem,'facealpha',0.5,'linestyle','none'); % plot the mesh

hold on;
% cut the mesh at value=-4
[cutpos,cutvalue,facedata]=qmeshcut(elem(:,1:4),node(:,1:3),log10(abs(g12)),-4); 
patch('Vertices',cutpos,'Faces',facedata,'FaceVertexCData',cutvalue,'FaceColor','interp');

% cut the mesh at value=-4.5
[cutpos,cutvalue,facedata]=qmeshcut(elem(:,1:4),node(:,1:3),log10(abs(g12)),-4.5);
patch('Vertices',cutpos,'Faces',facedata,'FaceVertexCData',cutvalue,'FaceColor','interp');

