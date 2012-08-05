function [gridOUTPUT,varargout] = VOXELISE(gridX,gridY,gridZ,varargin)
% VOXELISE  Voxelise a 3D triangular-polygon mesh
%==========================================================================
% FILENAME:          VOXELISE.m
% AUTHOR:            Adam H. Aitkenhead
% INSTITUTION:       The Christie NHS Foundation Trust
% CONTACT:           adam.aitkenhead@physics.cr.man.ac.uk
% DATE:              10th May 2010
% PURPOSE:           Voxelise a 3D triangular-polygon mesh.  The mesh must
%                    be closed (ie. watertight).
%
% USAGE:             [gridOUTPUT,gridCOx,gridCOy,gridCOz] = VOXELISE(gridX,gridY,gridZ,STLin,raydirection)
%             or...  [gridOUTPUT,gridCOx,gridCOy,gridCOz] = VOXELISE(gridX,gridY,gridZ,meshFV,raydirection)
%             or...  [gridOUTPUT,gridCOx,gridCOy,gridCOz] = VOXELISE(gridX,gridY,gridZ,meshX,meshY,meshZ,raydirection)
%             or...  [gridOUTPUT,gridCOx,gridCOy,gridCOz] = VOXELISE(gridX,gridY,gridZ,meshXYZ,raydirection)
%
% INPUT PARAMETERS:
%
%     gridX   - Mandatory - 1xP array     - List of the grid X coordinates. 
%                           OR an integer - Number of voxels in the grid in the X direction.
%
%     gridY   - Mandatory - 1xQ array     - List of the grid Y coordinates.
%                           OR an integer - Number of voxels in the grid in the Y direction.
%
%     gridZ   - Mandatory - 1xR array     - List of the grid Z coordinates.
%                           OR an integer - Number of voxels in the grid in the Z direction.
%
%     STLin   - Optional - string      - Filename of the STL file.
%
%     meshFV  - Optional - structure   - Structure containing the faces and
%                                        vertices of the mesh, in the same
%                                        format as that produced by the
%                                        isosurface command.  
%
%     meshX   - Optional - 3xN array   - List of the mesh X coordinates for the 3 vertices of each of the N triangular patches
%     meshY   - Optional - 3xN array   - List of the mesh Y coordinates for the 3 vertices of each of the N triangular patches
%     meshZ   - Optional - 3xN array   - List of the mesh Z coordinates for the 3 vertices of each of the N triangular patches
%
%     meshXYZ - Optional - Nx3x3 array - The vertex coordinates for each facet, with:
%                                        1 row for each facet
%                                        3 columns for the x,y,z coordinates
%                                        3 pages for the three vertices
%
%     raydirection - Optional - String - Defines the directions in which
%                                        ray-tracing is performed.  The
%                                        default is 'xyz', which traces in
%                                        the x,y,z directions and combines
%                                        the results.
%
% OUTPUT PARAMETERS:
%
%     gridOUTPUT - Mandatory - PxQxR logical array - Voxelised data (1=>Inside the mesh, 0=>Outside the mesh)
%
%     gridCOx    - Optional - 1xP array - List of the grid X coordinates.
%     gridCOy    - Optional - 1xQ array - List of the grid Y coordinates.
%     gridCOz    - Optional - 1xR array - List of the grid Z coordinates.
%
% EXAMPLES:
%
%     To voxelise an STL file:
%     >>  [gridOUTPUT] = VOXELISE(gridX,gridY,gridZ,STLin)
% 
%     To voxelise a mesh defined by a structure containing the faces and vertices:
%     >>  [gridOUTPUT] = VOXELISE(gridX,gridY,gridZ,meshFV)
% 
%     To voxelise a mesh where the x,y,z coordinates are defined by three 3xN arrays:
%     >>  [gridOUTPUT] = VOXELISE(gridX,gridY,gridZ,meshX,meshY,meshZ)
%
%     To voxelise a mesh defined by a single Nx3x3 array:
%     >>  [gridOUTPUT] = VOXELISE(gridX,gridY,gridZ,meshXYZ)
%
%     To also output the lists of X,Y,Z coordinates:
%     >>  [gridOUTPUT,gridCOx,gridCOy,gridCOz] = VOXELISE(gridX,gridY,gridZ,STLin)
%
%     To use ray-tracing in only the z-direction:
%     >>  [gridOUTPUT] = VOXELISE(gridX,gridY,gridZ,STLin,'z')
%
% NOTES:
%
%   - Defining raydirection='xyz' means that the mesh is ray-traced in each
%     of the x,y,z directions, with the overall result being a combination
%     of the result from each direction.  This gives the most reliable
%     result at the expense of computation time.
%   - Tracing in only one direction (eg. raydirection='z') is faster, but
%     can potentially lead to minor errors where a ray exactly crosses
%     several facet edges.
%
% REFERENCES:
%
%   - This code uses a ray intersection method similar to that described by:
%     Patil S and Ravi B.  Voxel-based representation, display and
%     thickness analysis of intricate shapes. Ninth International
%     Conference on Computer Aided Design and Computer Graphics (CAD/CG
%     2005)
%==========================================================================

%==========================================================================
% VERSION  USER  CHANGES
% -------  ----  -------
% 100510   AHA   Original version.
% 100514   AHA   Now works with non-STL input.  Changes also provide a
%                significant speed improvement.
% 100520   AHA   Now optionally output the grid x,y,z coordinates.
%                Robustness also improved.
% 100614   AHA   Define gridOUTPUT as a logical array to improve memory
%                efficiency.
% 100615   AHA   Reworked the ray interpolation code to correctly handle
%                logical arrays.
% 100623   AHA   Enable ray-tracing in any combination of the x,y,z
%                directions.
% 100628   AHA   Allow input to be a structure containing [faces,vertices]
%                data, similar to the type of structure output by
%                isosurface.
% 100907   AHA   Now allow the grid to be smaller than the mesh dimensions.
% 101126   AHA   Simplified code, slight speed improvement, more robust.
%                Changed handling of automatic grid generation to reduce
%                chance of artefacts.
% 101201   AHA   Fixed bug in automatic grid generation.
% 110303   AHA   Improved method of finding which mesh facets can possibly
%                be crossed by each ray.  Up to 80% reduction in run-time.
%==========================================================================


%======================================================
% CHECK THE REQUIRED NUMBER OF OUTPUT PARAMETERS
%======================================================

if nargout~=1 && nargout~=4
  error('Incorrect number of output arguments.')
end

%======================================================
% READ INPUT PARAMETERS
%======================================================

% Read the ray direction if defined by the user, and remove this from the
% list of input arguments.  This makes it to make it easier to extract the
% mesh data from the input arguments in the subsequent step.

if ischar(varargin{end}) && max(strcmpi(varargin{end},{'x','y','z','xy','xz','yx','yz','zx','zy','xyz','xzy','yxz','yzx','zxy','zyx'}))
  raydirection    = lower(varargin{end});
  varargin        = varargin(1:nargin-4);
  narginremaining = nargin-1;
else
  raydirection    = 'xyz';    %Set the default value if none is defined by the user
  narginremaining = nargin;
end

% Whatever the input mesh format is, it is converted to an Nx3x3 array
% defining the vertex positions for each facet, with 1 row for each facet,
% 3 cols for the x,y,z coordinates, and 3 pages for the three vertices.

if narginremaining==4
  
  if isstruct(varargin{1})==1
     
    faces  = varargin{1}.faces;
    vertex = varargin{1}.vertices;
    
    meshXYZ = zeros(size(faces,1),3,3);
    for loopa = 1:size(faces,1)
      meshXYZ(loopa,:,1) = vertex(faces(loopa,1),:);
      meshXYZ(loopa,:,2) = vertex(faces(loopa,2),:);
      meshXYZ(loopa,:,3) = vertex(faces(loopa,3),:);
    end
  
  elseif ischar(varargin{1})
    meshXYZ = READ_stl(varargin{1});
    
  else
    meshXYZ = varargin{1};
    
  end

elseif narginremaining==6

  meshX = varargin{1};
  meshY = varargin{2};
  meshZ = varargin{3};
 
  meshXYZ = zeros( size(meshX,2) , 3 , size(meshX,1) );

  meshXYZ(:,1,:) = reshape(meshX',size(meshX,2),1,3);
  meshXYZ(:,2,:) = reshape(meshY',size(meshY,2),1,3);
  meshXYZ(:,3,:) = reshape(meshZ',size(meshZ,2),1,3);
  
else
  
  error('Incorrect number of input arguments.')
  
end

%======================================================
% IDENTIFY THE MIN AND MAX X,Y,Z COORDINATES OF THE POLYGON MESH
%======================================================

meshXmin = min(min(meshXYZ(:,1,:)));
meshXmax = max(max(meshXYZ(:,1,:)));
meshYmin = min(min(meshXYZ(:,2,:)));
meshYmax = max(max(meshXYZ(:,2,:)));
meshZmin = min(min(meshXYZ(:,3,:)));
meshZmax = max(max(meshXYZ(:,3,:)));

%======================================================
% CHECK THE DIMENSIONS OF THE 3D OUTPUT GRID
%======================================================
% The output grid will be defined by the coordinates in gridCOx, gridCOy, gridCOz

if numel(gridX)>1
  if size(gridX,1)>size(gridX,2)   %gridX should be a row vector rather than a column vector
    gridCOx = gridX';
  else
    gridCOx = gridX;
  end
elseif numel(gridX)==1 && gridX==1   %If gridX is a single integer (rather than a vector) and is equal to 1
  gridCOx   = (meshXmin+meshXmax)/2;
elseif numel(gridX)==1 && rem(gridX,1)==0   %If gridX is a single integer (rather than a vector) then automatically create the list of x coordinates
  voxwidth  = (meshXmax-meshXmin)/(gridX+1/2);
  gridCOx   = [meshXmin+voxwidth/2:voxwidth:meshXmax-voxwidth/2];
end

if numel(gridY)>1
  if size(gridY,1)>size(gridY,2)   %gridY should be a row vector rather than a column vector
    gridCOy = gridY';
  else
    gridCOy = gridY;
  end
elseif numel(gridY)==1 && gridY==1   %If gridX is a single integer (rather than a vector) and is equal to 1
  gridCOy   = (meshYmin+meshYmax)/2;
elseif numel(gridY)==1 && rem(gridY,1)==0   %If gridX is a single integer (rather than a vector) then automatically create the list of y coordinates
  voxwidth  = (meshYmax-meshYmin)/(gridY+1/2);
  gridCOy   = [meshYmin+voxwidth/2:voxwidth:meshYmax-voxwidth/2];
end

if numel(gridZ)>1
  if size(gridZ,1)>size(gridZ,2)   %gridZ should be a row vector rather than a column vector
    gridCOz = gridZ';
  else
    gridCOz = gridZ;
  end
elseif numel(gridZ)==1 && gridZ==1   %If gridX is a single integer (rather than a vector) and is equal to 1
  gridCOz   = (meshZmin+meshZmax)/2;
elseif numel(gridZ)==1 && rem(gridZ,1)==0   %If gridZ is a single integer (rather than a vector) then automatically create the list of z coordinates
  voxwidth  = (meshZmax-meshZmin)/(gridZ+1/2);
  gridCOz   = [meshZmin+voxwidth/2:voxwidth:meshZmax-voxwidth/2];
end

%Check that the output grid is large enough to cover the mesh:
if ~isempty(strfind(raydirection,'x'))  &&  (min(gridCOx)>meshXmin || max(gridCOx)<meshXmax)
  gridcheckX = 0;
  if min(gridCOx)>meshXmin
    gridCOx = [meshXmin,gridCOx];
    gridcheckX = gridcheckX+1;
  end
  if max(gridCOx)<meshXmax
    gridCOx = [gridCOx,meshXmax];
    gridcheckX = gridcheckX+2;
  end
elseif ~isempty(strfind(raydirection,'y'))  &&  (min(gridCOy)>meshYmin || max(gridCOy)<meshYmax)
  gridcheckY = 0;
  if min(gridCOy)>meshYmin
    gridCOy = [meshYmin,gridCOy];
    gridcheckY = gridcheckY+1;
  end
  if max(gridCOy)<meshYmax
    gridCOy = [gridCOy,meshYmax];
    gridcheckY = gridcheckY+2;
  end
elseif ~isempty(strfind(raydirection,'z'))  &&  (min(gridCOz)>meshZmin || max(gridCOz)<meshZmax)
  gridcheckZ = 0;
  if min(gridCOz)>meshZmin
    gridCOz = [meshZmin,gridCOz];
    gridcheckZ = gridcheckZ+1;
  end
  if max(gridCOz)<meshZmax
    gridCOz = [gridCOz,meshZmax];
    gridcheckZ = gridcheckZ+2;
  end
end

%======================================================
% VOXELISE USING THE USER DEFINED RAY DIRECTION(S)
%======================================================

%Count the number of voxels in each direction:
voxcountX = numel(gridCOx);
voxcountY = numel(gridCOy);
voxcountZ = numel(gridCOz);

% Prepare logical array to hold the voxelised data:
gridOUTPUT = false( voxcountX,voxcountY,voxcountZ,numel(raydirection) );
countdirections = 0;

if strfind(raydirection,'x')
  countdirections = countdirections + 1;
  gridOUTPUT(:,:,:,countdirections) = permute( VOXELISEinternal(gridCOy,gridCOz,gridCOx,meshXYZ(:,[2,3,1],:)) ,[3,1,2] );
end

if strfind(raydirection,'y')
  countdirections = countdirections + 1;
  gridOUTPUT(:,:,:,countdirections) = permute( VOXELISEinternal(gridCOz,gridCOx,gridCOy,meshXYZ(:,[3,1,2],:)) ,[2,3,1] );
end

if strfind(raydirection,'z')
  countdirections = countdirections + 1;
  gridOUTPUT(:,:,:,countdirections) = VOXELISEinternal(gridCOx,gridCOy,gridCOz,meshXYZ);
end

% Combine the results of each ray-tracing direction:
if numel(raydirection)>1
  gridOUTPUT = sum(gridOUTPUT,4)>=numel(raydirection)/2;
end

%======================================================
% RETURN THE OUTPUT GRID TO THE SIZE REQUIRED BY THE USER (IF IT WAS CHANGED EARLIER)
%======================================================

if exist('gridcheckX','var')
  if gridcheckX == 1
    gridOUTPUT = gridOUTPUT(2:end,:,:);
    gridCOx    = gridCOx(2:end);
  elseif gridcheckX == 2
    gridOUTPUT = gridOUTPUT(1:end-1,:,:);
    gridCOx    = gridCOx(1:end-1);
  elseif gridcheckX == 3
    gridOUTPUT = gridOUTPUT(2:end-1,:,:);
    gridCOx    = gridCOx(2:end-1);
  end
end
if exist('gridcheckY','var')
  if gridcheckY == 1
    gridOUTPUT = gridOUTPUT(:,2:end,:);
    gridCOy    = gridCOy(2:end);
  elseif gridcheckY == 2
    gridOUTPUT = gridOUTPUT(:,1:end-1,:);
    gridCOy    = gridCOy(1:end-1);
  elseif gridcheckY == 3
    gridOUTPUT = gridOUTPUT(:,2:end-1,:);
    gridCOy    = gridCOy(2:end-1);
  end
end
if exist('gridcheckZ','var')
  if gridcheckZ == 1
    gridOUTPUT = gridOUTPUT(:,:,2:end);
    gridCOz    = gridCOz(2:end);
  elseif gridcheckZ == 2
    gridOUTPUT = gridOUTPUT(:,:,1:end-1);
    gridCOz    = gridCOz(1:end-1);
  elseif gridcheckZ == 3
    gridOUTPUT = gridOUTPUT(:,:,2:end-1);
    gridCOz    = gridCOz(2:end-1);
  end
end

%======================================================
% PREPARE THE OUTPUT PARAMETERS
%======================================================

if nargout==4
  varargout(1) = {gridCOx};
  varargout(2) = {gridCOy};
  varargout(3) = {gridCOz};
end

end %function
%==========================================================================





%==========================================================================
function [gridOUTPUT] = VOXELISEinternal(gridCOx,gridCOy,gridCOz,meshXYZ)

%Count the number of voxels in each direction:
voxcountX = numel(gridCOx);
voxcountY = numel(gridCOy);
voxcountZ = numel(gridCOz);

% Prepare logical array to hold the voxelised data:
gridOUTPUT = false(voxcountX,voxcountY,voxcountZ);

%Identify the min and max x,y coordinates (cm) of the mesh:
meshXmin = min(min(meshXYZ(:,1,:)));
meshXmax = max(max(meshXYZ(:,1,:)));
meshYmin = min(min(meshXYZ(:,2,:)));
meshYmax = max(max(meshXYZ(:,2,:)));
meshZmin = min(min(meshXYZ(:,3,:)));
meshZmax = max(max(meshXYZ(:,3,:)));

%Identify the min and max x,y coordinates (pixels) of the mesh:
meshXminp = find(abs(gridCOx-meshXmin)==min(abs(gridCOx-meshXmin)));
meshXmaxp = find(abs(gridCOx-meshXmax)==min(abs(gridCOx-meshXmax)));
meshYminp = find(abs(gridCOy-meshYmin)==min(abs(gridCOy-meshYmin)));
meshYmaxp = find(abs(gridCOy-meshYmax)==min(abs(gridCOy-meshYmax)));

%Make sure min < max for the mesh coordinates:
if meshXminp > meshXmaxp
  [meshXminp,meshXmaxp] = deal(meshXmaxp,meshXminp);
end %if
if meshYminp > meshYmaxp
  [meshYminp,meshYmaxp] = deal(meshYmaxp,meshYminp);
end %if

%Identify the min and max x,y,z coordinates of each facet:
meshXYZmin = min(meshXYZ,[],3);
meshXYZmax = max(meshXYZ,[],3);

%======================================================
% TURN OFF DIVIDE-BY-ZERO WARNINGS
%======================================================
%This prevents the Y1predicted, Y2predicted, Y3predicted and YRpredicted
%calculations creating divide-by-zero warnings.  Suppressing these warnings
%doesn't affect the code, because only the sign of the result is important.
%That is, 'Inf' and '-Inf' results are ok.
%The warning will be returned to its original state at the end of the code.
warningrestorestate = warning('query', 'MATLAB:divideByZero');
warning off MATLAB:divideByZero

%======================================================
% VOXELISE THE MESH
%======================================================

correctionLIST = [];   %Prepare to record all rays that fail the voxelisation.  This array is built on-the-fly, but since
                       %it ought to be relatively small should not incur too much of a speed penalty.
                       
% Loop through each x,y pixel.
% The mesh will be voxelised by passing rays in the z-direction through
% each x,y pixel, and finding the locations where the rays cross the mesh.
for loopY = meshYminp:meshYmaxp
  
  % - 1a - Find which mesh facets could possibly be crossed by the ray:
  possibleCROSSLISTy = find( meshXYZmin(:,2)<gridCOy(loopY) & meshXYZmax(:,2)>gridCOy(loopY) );
  
  for loopX = meshXminp:meshXmaxp
      
    % - 1b - Find which mesh facets could possibly be crossed by the ray:
    possibleCROSSLIST = possibleCROSSLISTy( meshXYZmin(possibleCROSSLISTy,1)<gridCOx(loopX) & meshXYZmax(possibleCROSSLISTy,1)>gridCOx(loopX) );

    if isempty(possibleCROSSLIST)==0  %Only continue the analysis if some nearby facets were actually identified
          
      % - 2 - For each facet, check if the ray really does cross the facet rather than just passing it close-by:
          
      % GENERAL METHOD:
      % 1. Take each edge of the facet in turn.
      % 2. Find the position of the opposing vertex to that edge.
      % 3. Find the position of the ray relative to that edge.
      % 4. Check if ray is on the same side of the edge as the opposing vertex.
      % 5. If this is true for all three edges, then the ray definitely passes through the facet.
      %
      % NOTES:
      % 1. If the ray crosses exactly on an edge, this is counted as crossing the facet.
      % 2. If a ray crosses exactly on a vertex, this is also taken into account.
      
      facetCROSSLIST = [];   %Prepare to record all facets which are crossed by the ray.  This array is built on-the-fly, but since
                             %it ought to be relatively small should not incur too much of a speed penalty.
                             
      for loopCHECKFACET = possibleCROSSLIST'
  
        %Check if ray crosses the facet.  This method is much (>>10 times) faster than using the built-in function 'inpolygon'.
        %Taking each edge of the facet in turn, check if the ray is on the same side as the opposing vertex.  If so, let testVn=1
        
        Y1predicted = meshXYZ(loopCHECKFACET,2,2) - ((meshXYZ(loopCHECKFACET,2,2)-meshXYZ(loopCHECKFACET,2,3)) * (meshXYZ(loopCHECKFACET,1,2)-meshXYZ(loopCHECKFACET,1,1))/(meshXYZ(loopCHECKFACET,1,2)-meshXYZ(loopCHECKFACET,1,3)));
        YRpredicted = meshXYZ(loopCHECKFACET,2,2) - ((meshXYZ(loopCHECKFACET,2,2)-meshXYZ(loopCHECKFACET,2,3)) * (meshXYZ(loopCHECKFACET,1,2)-gridCOx(loopX))/(meshXYZ(loopCHECKFACET,1,2)-meshXYZ(loopCHECKFACET,1,3)));
        if (Y1predicted > meshXYZ(loopCHECKFACET,2,1) && YRpredicted > gridCOy(loopY)) || (Y1predicted < meshXYZ(loopCHECKFACET,2,1) && YRpredicted < gridCOy(loopY))
          testV1 = 1;   %The ray is on the same side of the 2-3 edge as the 1st vertex.
        else
          testV1 = 0;   %The ray is on the opposite side of the 2-3 edge to the 1st vertex.
        end %if
          
        Y2predicted = meshXYZ(loopCHECKFACET,2,3) - ((meshXYZ(loopCHECKFACET,2,3)-meshXYZ(loopCHECKFACET,2,1)) * (meshXYZ(loopCHECKFACET,1,3)-meshXYZ(loopCHECKFACET,1,2))/(meshXYZ(loopCHECKFACET,1,3)-meshXYZ(loopCHECKFACET,1,1)));
        YRpredicted = meshXYZ(loopCHECKFACET,2,3) - ((meshXYZ(loopCHECKFACET,2,3)-meshXYZ(loopCHECKFACET,2,1)) * (meshXYZ(loopCHECKFACET,1,3)-gridCOx(loopX))/(meshXYZ(loopCHECKFACET,1,3)-meshXYZ(loopCHECKFACET,1,1)));
        if (Y2predicted > meshXYZ(loopCHECKFACET,2,2) && YRpredicted > gridCOy(loopY)) || (Y2predicted < meshXYZ(loopCHECKFACET,2,2) && YRpredicted < gridCOy(loopY))
          testV2 = 1;   %The ray is on the same side of the 3-1 edge as the 2nd vertex.
        else
          testV2 = 0;   %The ray is on the opposite side of the 3-1 edge to the 2nd vertex.
        end %if

        Y3predicted = meshXYZ(loopCHECKFACET,2,1) - ((meshXYZ(loopCHECKFACET,2,1)-meshXYZ(loopCHECKFACET,2,2)) * (meshXYZ(loopCHECKFACET,1,1)-meshXYZ(loopCHECKFACET,1,3))/(meshXYZ(loopCHECKFACET,1,1)-meshXYZ(loopCHECKFACET,1,2)));
        YRpredicted = meshXYZ(loopCHECKFACET,2,1) - ((meshXYZ(loopCHECKFACET,2,1)-meshXYZ(loopCHECKFACET,2,2)) * (meshXYZ(loopCHECKFACET,1,1)-gridCOx(loopX))/(meshXYZ(loopCHECKFACET,1,1)-meshXYZ(loopCHECKFACET,1,2)));
        if (Y3predicted > meshXYZ(loopCHECKFACET,2,3) && YRpredicted > gridCOy(loopY)) || (Y3predicted < meshXYZ(loopCHECKFACET,2,3) && YRpredicted < gridCOy(loopY))
          testV3 = 1;   %The ray is on the same side of the 1-2 edge as the 3rd vertex.
        else
          testV3 = 0;   %The ray is on the opposite side of the 1-2 edge to the 3rd vertex.
        end %if
        
        %The ray passes through the facet if it is on the correct side of
        %all 3 edges OR it crosses exactly one edge.
        if testV1==1 && testV2==1 && testV3==1
          facetCROSSLIST = [facetCROSSLIST,loopCHECKFACET];
        end %if
      
      end %for
    

      % - 3 - Find the z coordinate of the locations where the ray crosses each facet:

      gridCOzCROSS = zeros(size(facetCROSSLIST));
      for loopFINDZ = facetCROSSLIST

        % METHOD:
        % 1. Define the equation describing the plane of the facet.  For a
        % more detailed outline of the maths, see:
        % http://local.wasp.uwa.edu.au/~pbourke/geometry/planeeq/
        %    Ax + By + Cz + D = 0
        %    where  A = y1 (z2 - z3) + y2 (z3 - z1) + y3 (z1 - z2)
        %           B = z1 (x2 - x3) + z2 (x3 - x1) + z3 (x1 - x2)
        %           C = x1 (y2 - y3) + x2 (y3 - y1) + x3 (y1 - y2)
        %           D = - x1 (y2 z3 - y3 z2) - x2 (y3 z1 - y1 z3) - x3 (y1 z2 - y2 z1)
        % 2. For the x and y coordinates of the ray, solve these equations to find the z coordinate in this plane.

        planecoA = meshXYZ(loopFINDZ,2,1)*(meshXYZ(loopFINDZ,3,2)-meshXYZ(loopFINDZ,3,3)) + meshXYZ(loopFINDZ,2,2)*(meshXYZ(loopFINDZ,3,3)-meshXYZ(loopFINDZ,3,1)) + meshXYZ(loopFINDZ,2,3)*(meshXYZ(loopFINDZ,3,1)-meshXYZ(loopFINDZ,3,2));
        planecoB = meshXYZ(loopFINDZ,3,1)*(meshXYZ(loopFINDZ,1,2)-meshXYZ(loopFINDZ,1,3)) + meshXYZ(loopFINDZ,3,2)*(meshXYZ(loopFINDZ,1,3)-meshXYZ(loopFINDZ,1,1)) + meshXYZ(loopFINDZ,3,3)*(meshXYZ(loopFINDZ,1,1)-meshXYZ(loopFINDZ,1,2)); 
        planecoC = meshXYZ(loopFINDZ,1,1)*(meshXYZ(loopFINDZ,2,2)-meshXYZ(loopFINDZ,2,3)) + meshXYZ(loopFINDZ,1,2)*(meshXYZ(loopFINDZ,2,3)-meshXYZ(loopFINDZ,2,1)) + meshXYZ(loopFINDZ,1,3)*(meshXYZ(loopFINDZ,2,1)-meshXYZ(loopFINDZ,2,2));
        planecoD = - meshXYZ(loopFINDZ,1,1)*(meshXYZ(loopFINDZ,2,2)*meshXYZ(loopFINDZ,3,3)-meshXYZ(loopFINDZ,2,3)*meshXYZ(loopFINDZ,3,2)) - meshXYZ(loopFINDZ,1,2)*(meshXYZ(loopFINDZ,2,3)*meshXYZ(loopFINDZ,3,1)-meshXYZ(loopFINDZ,2,1)*meshXYZ(loopFINDZ,3,3)) - meshXYZ(loopFINDZ,1,3)*(meshXYZ(loopFINDZ,2,1)*meshXYZ(loopFINDZ,3,2)-meshXYZ(loopFINDZ,2,2)*meshXYZ(loopFINDZ,3,1));

        if abs(planecoC) < 1e-14
          planecoC=0;
        end
      
        gridCOzCROSS(facetCROSSLIST==loopFINDZ) = (- planecoD - planecoA*gridCOx(loopX) - planecoB*gridCOy(loopY)) / planecoC;
        
      end %for

      %Remove values of gridCOzCROSS which are outside of the mesh limits (including a 1e-12 margin for error).
      gridCOzCROSS = gridCOzCROSS( gridCOzCROSS>=meshZmin-1e-12 & gridCOzCROSS<=meshZmax+1e-12 );
      
      %Round gridCOzCROSS to remove any rounding errors, and take only the unique values:
      gridCOzCROSS = round(gridCOzCROSS*1e12)/1e12;
      gridCOzCROSS = unique(gridCOzCROSS);
      
      
      % - 4 - Label as being inside the mesh all the voxels that the ray passes through after crossing one facet before crossing another facet:

      if rem(numel(gridCOzCROSS),2)==0  % Only rays which cross an even number of facets are voxelised

        for loopASSIGN = 1:(numel(gridCOzCROSS)/2)
          voxelsINSIDE = (gridCOz>gridCOzCROSS(2*loopASSIGN-1) & gridCOz<gridCOzCROSS(2*loopASSIGN));
          gridOUTPUT(loopX,loopY,voxelsINSIDE) = 1;
        end %for
        
      elseif numel(gridCOzCROSS)~=0    % Remaining rays which meet the mesh in some way are not voxelised, but are labelled for correction later.
          
        correctionLIST = [ correctionLIST; loopX,loopY ];
        
      end %if

    end %if

  end %for
end %for

%======================================================
% RESTORE DIVIDE-BY-ZERO WARNINGS TO THE ORIGINAL STATE
%======================================================

warning(warningrestorestate)

%======================================================
% USE INTERPOLATION TO FILL IN THE RAYS WHICH COULD NOT BE VOXELISED
%======================================================
%For rays where the voxelisation did not give a clear result, the ray is
%computed by interpolating from the surrounding rays.
countCORRECTIONLIST = size(correctionLIST,1);

if countCORRECTIONLIST>0
    
  %If necessary, add a one-pixel border around the x and y edges of the
  %array.  This prevents an error if the code tries to interpolate a ray at
  %the edge of the x,y grid.
  if min(correctionLIST(:,1))==1 || max(correctionLIST(:,1))==numel(gridCOx) || min(correctionLIST(:,2))==1 || max(correctionLIST(:,2))==numel(gridCOy)
    gridOUTPUT     = [zeros(1,voxcountY+2,voxcountZ);zeros(voxcountX,1,voxcountZ),gridOUTPUT,zeros(voxcountX,1,voxcountZ);zeros(1,voxcountY+2,voxcountZ)];
    correctionLIST = correctionLIST + 1;
  end
  
  for loopC = 1:countCORRECTIONLIST
    voxelsforcorrection = squeeze( sum( [ gridOUTPUT(correctionLIST(loopC,1)-1,correctionLIST(loopC,2)-1,:) ,...
                                          gridOUTPUT(correctionLIST(loopC,1)-1,correctionLIST(loopC,2),:)   ,...
                                          gridOUTPUT(correctionLIST(loopC,1)-1,correctionLIST(loopC,2)+1,:) ,...
                                          gridOUTPUT(correctionLIST(loopC,1),correctionLIST(loopC,2)-1,:)   ,...
                                          gridOUTPUT(correctionLIST(loopC,1),correctionLIST(loopC,2)+1,:)   ,...
                                          gridOUTPUT(correctionLIST(loopC,1)+1,correctionLIST(loopC,2)-1,:) ,...
                                          gridOUTPUT(correctionLIST(loopC,1)+1,correctionLIST(loopC,2),:)   ,...
                                          gridOUTPUT(correctionLIST(loopC,1)+1,correctionLIST(loopC,2)+1,:) ,...
                                         ] ) );
    voxelsforcorrection = (voxelsforcorrection>=4);
    gridOUTPUT(correctionLIST(loopC,1),correctionLIST(loopC,2),voxelsforcorrection) = 1;
  end %for

  %Remove the one-pixel border surrounding the array, if this was added
  %previously.
  if size(gridOUTPUT,1)>numel(gridCOx) || size(gridOUTPUT,2)>numel(gridCOy)
    gridOUTPUT = gridOUTPUT(2:end-1,2:end-1,:);
  end
  
end %if

%disp([' Ray tracing result: ',num2str(countCORRECTIONLIST),' rays (',num2str(countCORRECTIONLIST/(voxcountX*voxcountY)*100,'%5.1f'),'% of all rays) exactly crossed a facet edge and had to be computed by interpolation.'])

end %function
%==========================================================================
