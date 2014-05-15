% test_itk_tri_rasterization.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cube mesh
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cube mesh vertices, a bit away from voxel centres to avoid inconsistent
% results (see next section)
x = [
    0.5 -0.1 0
    0.5 -0.1 1.01
    0.5 1.1 0
    0.5 1.1 1.01
    1.505 -0.1 0
    1.505 -0.1 1.01
    1.505 1.1 0
    1.505 1.1 1.01
    ];


% apply alphashape
[~, as] = alphavol(x, 1.3);
tri = as.bnd;
clear as

% plot mesh
hold off
trimesh(tri, x(:,1), x(:,2), x(:,3))
axis([0.5 1.5 0 1 0 1])

% compute parameters of output image
res = [.2 .1 .1]; % r, c, s
origin = [0 -.5 -.5]; % x, y, z
final = [2 1.5 1.5]; % x, y, z
sz = 1 + round((final([2 1 3]) - origin([2 1 3]))./res); % r, c, s

% centers of voxels in the image
cx = linspace(origin(1), final(1), sz(2));
cy = linspace(origin(2), final(2), sz(1));
cz = linspace(origin(3), final(3), sz(3));

% rasterize the mesh
im = itk_tri_rasterization(tri, x, res, sz, origin);

% plot segmentation mask, and overlay mesh
subplot(1, 2, 1)
hold off
imagesc([origin(1) final(1)], [origin(2) final(2)], im(:, :, 6))
hold on
for I = 1:length(cx)
    plot(cx(I)*[1 1], [origin(2) final(2)])
end
for I = 1:length(cy)
    plot([origin(1) final(1)], cy(I)*[1 1])
end
xmin = min(x(:, 1));
xmax = max(x(:, 1));
ymin = min(x(:, 2));
ymax = max(x(:, 2));
plot([xmin xmax xmax xmin xmin], [ymin ymin ymax ymax ymin], 'w')
xlabel('x')
ylabel('y')
axis xy

% plot a vertical cut
subplot(1, 2, 2)
hold off
% imagesc([origin(1) final(1)], [origin(3) final(3)], squeeze(im(:, 11, :)))
imagesc([origin(2) final(2)], [origin(3) final(3)], squeeze(im(:, 11, :))')
hold on
for I = 1:length(cy)
    plot(cy(I)*[1 1], [origin(3) final(3)])
end
for I = 1:length(cz)
    plot([origin(2) final(2)], cz(I)*[1 1])
end
zmin = min(x(:, 3));
zmax = max(x(:, 3));
plot([ymin ymax ymax ymin ymin], [zmin zmin zmax zmax zmin], 'w')
xlabel('y')
ylabel('z')
axis xy

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Note that voxelize offers inconsistent results if the mesh exactly goes
%% through a voxel centre
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% depending on whether the voxel is at the bottom, top, left or right of
% the segmentation, it will be or not set to 1. This can be seen using the
% following set of vertices.

% create a surface mesh that is a cube
x = [
    0.5 0 0
    0.5 0 1
    0.5 1 0
    0.5 1 1
    1.5 0 0
    1.5 0 1
    1.5 1 0
    1.5 1 1
    ];

