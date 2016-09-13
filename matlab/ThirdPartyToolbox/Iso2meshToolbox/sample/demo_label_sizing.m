%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% specify maximum cell sizes for different labels
% (this is only valid when using 'cgalmesh' with v2m or vol2mesh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create a dummy image

rad = 10:10:50;
maxRad = max(rad);
box = zeros(maxRad+1,maxRad+1,maxRad+1);
image = box;
box(1,1,1) = 1;
%DT = bwdist(box,'euclidean');
[ix,iy,iz]=meshgrid(1:maxRad+1,1:maxRad+1,1:maxRad+1);
DT = sqrt((ix - 1).*(ix -1) + (iy - 1).*(iy -1) + (iz - 1).*(iz -1));

for i = 1:size(rad,2)
    ball = DT<=rad(i);
    image(ball~=0) = image(ball~=0) + ball(ball~=0);
end

image = uint8(image);

% mesh the domain with different sizing options

figure;
maxvol='1';
[no,el]=v2m(image,[],5,maxvol,'cgalmesh');
subplot(221);
plotmesh(no(:,1:3),el,'x-y<0');
title('a single scalar sets cell size for all labels');

maxvol='1=2:2=1:3=2:4=1';
[no,el]=v2m(image,[],5,maxvol,'cgalmesh');
subplot(222);
plotmesh(no(:,1:3),el,'x-y<0');
title(sprintf('maxvol is "%s"',maxvol));

maxvol='2:1:2:1';
[no,el]=v2m(image,[],5,maxvol,'cgalmesh');
subplot(224);
plotmesh(no(:,1:3),el,'x-y<0');
title(sprintf('maxvol is "%s", same as above',maxvol));

maxvol='3=2:1:0.5';
[no,el]=v2m(image,[],5,maxvol,'cgalmesh');
subplot(223);
plotmesh(no(:,1:3),el,'x-y<0');
title(sprintf('maxvol is "%s", same to 3=2:4=1:5=0.5',maxvol));

