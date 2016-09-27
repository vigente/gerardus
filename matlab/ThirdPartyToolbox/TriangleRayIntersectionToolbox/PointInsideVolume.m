function inside = PointInsideVolume(point, faces, vertices)
%% Point within the volume test
% TriangleRayIntersection is a low level function which can be used to
% solve higher level problems. For example a test to see if point is inside
% or outside of a halume defined by a continous surface:

%% chack input
if nargin==2 && isa(faces, 'triangulation')
  [faces, vertices] = freeBoundary(faces);
end

switch size(point,2)
  case 2 % 2D case
    xv = vertices(faces', 1); % define polygon
    yv = vertices(faces', 2);
    inside = inpolygon(point(:,1), point(:,2), xv, yv);
  case 3
    eps    = 1e-5;
    vert1  = vertices(faces(:,1),:);
    vert2  = vertices(faces(:,2),:);
    vert3  = vertices(faces(:,3),:);
    inside = false(size(point,1),1);
    for iPoint = 1:size(point,1)
      certain = 0;
      while ~certain
        dir = rand(1,3)-0.5; % pick random direction
        [intersect, ~, u, v] = TriangleRayIntersection(point(iPoint,:), ...
          dir, vert1, vert2, vert3, 'border', 'inclusive');
        nIntersect = sum(intersect);    % number of intersections
        inside(iPoint) = mod(nIntersect,2)>0; % inside if odd number of intersections
        % make sure ray stays away fron surface triangle edges
        bary = [u, v, 1-u-v];
        bary = bary(intersect,:);
        certain = all( min(abs(bary),[], 2)>eps );
      end
    end
end