%% test untanglement with two free vertices

% 5 points in the perimeter, 2 free vertices, possible untangled solution
x = [
    0 3
    1 0
    3 1
    4 .5
    3.5 4
    1.5 2
    3 2
    ];

isFree = logical([0 0 0 0 0 1 1]');

tri = [
    6 1 2
    6 2 3
    7 3 4
    4 5 7
    5 1 7
    7 1 6
    7 6 3
    ];

% plot mesh
d = dmatrix_mesh(tri);
hold off
gplot(d, x)
hold on
plot(x(:, 1), x(:, 2), 'o')

% compute untangled coordinates for the free vertices
u = vertices_untangle(tri, x, isFree);

% check the solution
y = x;
y(isFree, :) = u;

% plot the solution
hold on
gplot(d, y, 'r')
plot(y(:, 1), y(:, 2), 'or')

%% test untanglement with three free vertices

% 5 points in the perimeter, 3 free vertices, possible untangled solution
x = [
    0 3
    1 0
    3 1
    4 .5
    3.5 4
    1.5 2
    3 2
    2 1.25
    ];

isFree = logical([0 0 0 0 0 1 1 1]');

tri = [
    6 1 2
    8 6 2
    8 2 3
    7 3 4
    4 5 7
    5 1 7
    7 1 6
    7 6 8
    7 8 3
    ];

% plot mesh
d = dmatrix_mesh(tri);
hold off
gplot(d, x)
hold on
plot(x(:, 1), x(:, 2), 'o')

% compute untangled coordinates for the free vertices
[u, fval, exitflag, info] = vertices_untangle(tri, x, isFree);

% check the solution
y = x;
y(isFree, :) = u;

% plot the solution
hold on
gplot(d, y, 'r')
plot(y(:, 1), y(:, 2), 'or')
