function [varargout] = CONVERT_meshformat(varargin)
%CONVERT_meshformat  Convert mesh data from array to faces,vertices format or vice versa
%==========================================================================
% FILENAME:          CONVERT_meshformat.m
% AUTHOR:            Adam H. Aitkenhead
% DATE:              17th August 2010
% PURPOSE:           Convert mesh data from array to faces,vertices format
%
% USAGE:             [faces,vertices] = CONVERT_meshformat(coordVERTICES)
%          ...or...  [coordVERTICES]  = CONVERT_meshformat(faces,vertices)
%==========================================================================

%==========================================================================
% VERSION  USER  CHANGES
% -------  ----  -------
% 100817   AHA   Original version
%==========================================================================


if nargin==2 && nargout==1

    
  faces  = varargin{1};
  vertex = varargin{2};
   
  meshXYZ = zeros(size(faces,1),3,3);
  for loopa = 1:size(faces,1)
    meshXYZ(loopa,:,1) = vertex(faces(loopa,1),:);
    meshXYZ(loopa,:,2) = vertex(faces(loopa,2),:);
    meshXYZ(loopa,:,3) = vertex(faces(loopa,3),:);
  end

  varargout(1) = {meshXYZ};
  
  
elseif nargin==1 && nargout==2

    
  meshXYZ = varargin{1};
  
  vertices = [meshXYZ(:,:,1);meshXYZ(:,:,2);meshXYZ(:,:,3)];
  vertices = unique(vertices,'rows');

  faces = zeros(size(meshXYZ,1),3);

  for loopF = 1:size(meshXYZ,1)
    for loopV = 1:3
        
      %[C,IA,vertref] = intersect(meshXYZ(loopF,:,loopV),vertices,'rows');
      %The following 3 lines are equivalent to the previous line, but are much faster:
      
      vertref = find(vertices(:,1)==meshXYZ(loopF,1,loopV));
      vertref = vertref(find(vertices(vertref,2)==meshXYZ(loopF,2,loopV)));
      vertref = vertref(find(vertices(vertref,3)==meshXYZ(loopF,3,loopV)));
      
      faces(loopF,loopV) = vertref;
    end
  end
  
  varargout(1) = {faces};
  varargout(2) = {vertices};
  
  
end
