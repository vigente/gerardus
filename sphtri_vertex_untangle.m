function latlon = sphtri_vertex_untangle(tri, latlon)


% DEBUG = false;
DEBUG = true;%%%%%%%%%%%%%%%

% number of vertices
Nv = size(latlon, 1);

% compute Cartesian coordinates of the sphere points
x = zeros(size(latlon, 1), 3);
[x(:, 1), x(:, 2), x(:, 3)] = sph2cart(latlon(:, 2), latlon(:, 1), 1);

% adjacency matrix for the spherical mesh
d = dmatrix_mesh(tri);

% find vertices that are tangled
nTangled = sphtri_vertex_istangled(tri, latlon);
nnz(nTangled)%%%%%%%%%%%%%%%%%%

%% one-by-one pre-untangling

% for vidx = find(nTangled)
%     
%     % get local neighbourhood and center it around lat=0, lon=0
%     [triloc, xloc, idxloc, R, idxtri] ...
%         = get_centered_local_neighbourhood(tri, x, vidx);
%     
%     % convert to latitude/longitude the centered local neighbourhood
%     [lonloc, latloc] = cart2sph(xloc(:, 1), xloc(:, 2), xloc(:, 3));
%     
%     % check whether this vertex is still tangled (it may have got
%     % untangled as a side-effect of another vertex being untangled). If
%     % it's not tangled anymore, we can just skip to the next one
%     if (~sphtri_vertex_istangled(triloc, [latloc, lonloc], idxloc))
%         continue
%     end
%     
%     % compute signed area of each triangle in the local neighbourhood
%     if (DEBUG)
%         %             vidx
%         %             a = trifacet_signed_area(triloc, [lonloc latloc])
%         
%         hold off
%         dloc = dmatrix_mesh(triloc);
%         gplot(dloc, [lonloc, latloc])
%         hold on
%         plot(0, 0, 'ro')
%     end
%     
%     % untangle the central vertex
%     [lonlatloc0, exitflag] ...
%         = vertex_untangle(triloc, [lonloc, latloc], idxloc);
%     lonloc(idxloc) = lonlatloc0(1);
%     latloc(idxloc) = lonlatloc0(2);
%     
%     % relocate central vertex only if linear programming found a solution
%     if (exitflag ~= 1)
%         continue
%     end
%     
%     if (DEBUG)
%         hold off
%         gplot(dloc, [lonloc, latloc])
%         hold on
%         plot(lonloc, latloc, 'o')
%         plot(lonloc(idxloc), latloc(idxloc), 'r*')
%     end
%     
%     % undo the centering rotation on the newly computed central vertex
%     xc = zeros(1, 3);
%     [xc(1), xc(2), xc(3)] ...
%         = sph2cart(lonloc(idxloc), latloc(idxloc), 1);
%     xc = xc * R';
%     [latlon(vidx, 2), latlon(vidx, 1)] = cart2sph(xc(1), xc(2), xc(3));
%     
%     % re-assess whether vertices are tangled
%     nTangled = sphtri_vertex_istangled(tri, latlon);
%     nnz(nTangled)%%%%%%%%%%%%%%%%%%
% 
% end

%% group untangling

% stop when there are no more tangled vertices
while (nnz(nTangled))
    
    % plot surface highlighting the tangled vertices
    if (DEBUG)
        hold off
        aux = zeros(size(latlon, 1), 3);
        [aux(:, 1), aux(:, 2), aux(:, 3)] ...
            = sph2cart(latlon(:, 2), latlon(:, 1), 1);
        trisurf(tri, aux(:, 1), aux(:, 2), aux(:, 3), nTangled)
        axis equal
    end
    
    % get the first tangled vertex from the list
    vidx = find(nTangled, 1, 'first');
    
    % find all connected vertices that are also tangled
    tangledboundaryidx = false(size(nTangled));
    tangledboundaryidx(vidx) = true;
    totalidx = tangledboundaryidx;
    while (nnz(tangledboundaryidx))
        
        % neighbours of the boundary tangled vertices
        nnidx = sum(d(tangledboundaryidx, :), 1) > 0;
        
%         find(nnidx)
        
        % remove the neighbours that are already in the local neighbourhood
        % (we want the boundary to expand outwards only)
        nnidx(totalidx) = false;
        
%         find(nnidx)
        
        % add these neighbours to the local neighbourhood
        totalidx(nnidx) = true;
        
%         find(totalidx)
        
        % vertices on the new boundary that are tangled
        tangledboundaryidx = nnidx & nTangled;
        
%         find(tangledboundaryidx)
        
    end
    
    % extract the local neighbourhood
    triloc = tri(sum(ismember(tri, find(totalidx)), 2) == 3, :);
    [triloc, xloc, iloc] = tri_squeeze(triloc, x);
    
    % count number of triangles each edge belongs to
    ed = sortrows(sort([triloc(:, 1:2); triloc(:, 2:3); triloc(:, [3 1])], 2));
    [edunique, ~, ied] = unique(ed, 'rows');
    count = hist(ied, 1:ied(end));
    
    % boundary edges are those that only belong to one triangle
    edboundary = edunique(count == 1, :);
    
    % boundary vertices
    vboundary = unique(edboundary);
%     
%     % get index of the central vertex in the local neighbourhood
%     [~, idxloc] = min((xloc(:, 1) - x(vidx, 1)).^2 ...
%         + (xloc(:, 2) - x(vidx, 2)).^2 ...
%         + (xloc(:, 3) - x(vidx, 3)).^2);
    
    % axis-angle representation of rotation that takes the current vertex
    % to [1, 0, 0]
    r = [cross(x(vidx, :), [1 0 0]) -acos(dot(x(vidx, :), [1 0 0]))];
    
    % 3D-rotation matrix representation of the axis-angle
    R = vrrotvec2mat(r);
    
    % rotate local neighbourhood so that the central voxel is on [1, 0, 0]
    xloc = xloc * R;
    
    % DEBUG
    dloc = dmatrix_mesh(triloc);
    [lonloc, latloc] = cart2sph(xloc(:, 1), xloc(:, 2), xloc(:, 3));
    hold off
    gplot(dloc, [lonloc, latloc])
    hold on
    plot(lonloc(vboundary), latloc(vboundary), 'ro')
    
    % list of free vertices
    isFree = true(length(lonloc), 1);
    isFree(vboundary) = false;
    
    % untangle vertices
    [u, fval, exitflag, info] = vertices_untangle(triloc, [lonloc, latloc], isFree);
    
    lonloc2 = lonloc;
    lonloc2(isFree) = u(:, 1);
    latloc2 = latloc;
    latloc2(isFree) = u(:, 2);

    % DEBUG
    hold off
    gplot(dloc, [lonloc2, latloc2])
    hold on
    plot(lonloc2(vboundary), latloc2(vboundary), 'ro')
    
    
    %%%%%%%%% UP To HERE
    
    
    % get a list of tangled vertices, sorted from less tangled to most
    % tangled
    nTangled(nTangled == 0) = Inf;
    [~, idxTangled] = sort(nTangled);
    idxTangled = idxTangled(~isinf(nTangled(idxTangled)));
    
    % try to untangle the vertices one by one
    for I = 1:length(idxTangled)
        
        % select vertex
        vidx = idxTangled(I);
        
    
end


% convert local neighbourhood to a mesh with each triangle in
% counter-clockwise orientation
% triloc = [
%           2 3 1
%           3 4 1
%           4 5 1
%           5 2 1
%          ]
triloc = [(2:length(vneigh)+1)' [3:length(vneigh)+1 2]' ones(length(vneigh), 1)];
xloc = [latlon(vidx, :) - latlonlocm;
    lonlatloc(aux(1:end-1), :)];

% compute new location for the tangled vertex such that is within the
% convex hull of neighbours
[latlon0, exitflag] = vertex_untangle(triloc, xloc);

% if the vertex had to be relocated and the linear programming method
% converged to a solution, we relocate it
if (exitflag==1)
    latlon(vidx, :) = latlon0 + latlonlocm;
    disp('relocated')
else
    disp('not relocated')
end

end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Auxiliary functions

function [triloc, xloc, idxloc, R, triidx] ...
    = get_centered_local_neighbourhood(tri, x, vidx)

% list of triangles that are adjacent to the current vertex
[triidx, ~] = find(ismember(tri, vidx));
    
% extract local neighbourhood from the whole mesh
[triloc, xloc] = tri_squeeze(tri(triidx, :), x);

% get index of the central vertex in the local neighbourhood
[~, idxloc] = min((xloc(:, 1) - x(vidx, 1)).^2 ...
    + (xloc(:, 2) - x(vidx, 2)).^2 ...
    + (xloc(:, 3) - x(vidx, 3)).^2);

% axis-angle representation of rotation that takes the current vertex
% to [1, 0, 0]
r = [cross(x(vidx, :), [1 0 0]) -acos(dot(x(vidx, :), [1 0 0]))];

% 3D-rotation matrix representation of the axis-angle
R = vrrotvec2mat(r);

% rotate local neighbourhood so that the central voxel is on [1, 0, 0]
xloc = xloc * R;

end

% % % function vn = get_sorted_neighbours(tri, d, v)
% % % 
% % % % unsorted neighbours of this vertex
% % % vn = find(d(:, v));
% % % 
% % % % put together the current vertex and its neighbours
% % % vall = [v; vn];
% % % 
% % % % list of triangles in the neighbourhood
% % % triloc = tri(sum(...
% % %     [ismember(tri(:, 1), vall), ...
% % %     ismember(tri(:, 2), vall), ...
% % %     ismember(tri(:, 3), vall)], 2)==3, : ...
% % %     );
% % % 
% % % % get ordering of the neighbours
% % % idx = sort_perim_vertices(d(vall, vall), 1);
% % % 
% % % % sort the neighbours
% % % vn = vall(idx);
% % % 
% % % end

% d is the sparse adjacency matrix with the connections between the
% perimeter vertices
%
% input v is a list of the vertices that we want to sort
%
% output v is the list of vertices sorted. Some of the input vertices may
% have been removed if they were producing loops in the boundary
% % % function v = sort_perim_vertices(d, v)
% % % 
% % % % keep only the connections between the boundary vertices
% % % V = size(d, 1);
% % % dloc = sparse(V, V);
% % % dloc(v, v) = d(v, v)>0;
% % % 
% % % % Dijkstra's shortest paths from one arbitrary vertex in the perimeter to
% % % % all the others
% % % [l, p] = dijkstra(dloc, v(1));
% % % 
% % % % arbitrarily choose one vertex at the antipodes in the perimeter 
% % % l(isinf(l)) = 0;
% % % [~, vb] = max(l);
% % % 
% % % % path from v(1) to vb
% % % branch1 = graphpred2path(p, vb);
% % % 
% % % % remove the connections in branch1, so that we can find the other branch
% % % dloc(sub2ind([V V], branch1(1:end-1)', branch1(2:end)')) = 0;
% % % dloc(sub2ind([V V], branch1(2:end)', branch1(1:end-1)')) = 0;
% % % 
% % % % recompute Dijkstra's shortest path from vb
% % % [~, p] = dijkstra(dloc, vb);
% % % 
% % % % 2nd path from va to vb
% % % branch2 = graphpred2path(p, v(1));
% % % 
% % % % combine both branches into final solution
% % % v = [branch1 branch2(end-1:-1:2)];
% % % 
% % % end
