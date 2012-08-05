function [coordVERTICES,varargout] = READ_stl(stlFILENAME,varargin)
% READ_stlascii  Read mesh data in the form of an <*.stl> file
%==========================================================================
% FILENAME:          READ_stl.m
% AUTHOR:            Adam H. Aitkenhead
% INSTITUTION:       The Christie NHS Foundation Trust
% CONTACT:           adam.aitkenhead@physics.cr.man.ac.uk
% DATE:              29th March 2010
% PURPOSE:           Read mesh data in the form of an <*.stl> file.
%
% USAGE:
%
%     [coordVERTICES,coordNORMALS,stlNAME] = READ_stl(stlFILENAME,stlFORMAT)
%
% INPUT PARAMETERS:
%
%     stlFILENAME   - String - Mandatory - The filename of the STL file.
%
%     stlFORMAT     - String - Optional  - The format of the STL file:
%                                        'ascii' or 'binary'
%
% OUTPUT PARAMETERS:
%
%     coordVERTICES - Nx3x3 array - Mandatory
%                                 - An array defining the vertex positions
%                                   for each of the N facets, with: 
%                                   1 row for each facet
%                                   3 cols for the x,y,z coordinates
%                                   3 pages for the three vertices
%
%     coordNORMALS  - Nx3 array   - Optional
%                                 - An array defining the normal vector for
%                                   each of the N facets, with: 
%                                   1 row for each facet
%                                   3 cols for the x,y,z components of the vector
%
%     stlNAME       - String      - Optional  - The name of the STL object.
%
%==========================================================================

%==========================================================================
% VERSION  USER  CHANGES
% -------  ----  -------
% 100329   AHA   Original version
% 100513   AHA   Totally reworked the code.  Now use textscan to read the
%                file all at once, rather than one line at a time with
%                fgetl.  Major speed improvment.
% 100623   AHA   Combined code which reads ascii STLS and code which reads
%                binary STLs into a single function.
% 101126   AHA   Small change to binary read code:  Now use fread instead
%                of fseek.  Gives large speed increase.
%==========================================================================

%==========================================================================
% STL ascii file format
%======================
% ASCII STL files have the following structure.  Technically each facet
% could be any 2D shape, but in practice only triangular facets tend to be
% used.  The present code ONLY works for meshes composed of triangular
% facets.
%
% solid object_name
% facet normal x y z
%   outer loop
%     vertex x y z
%     vertex x y z
%     vertex x y z
%   endloop
% endfacet
%
% <Repeat for all facets...>
%
% endsolid object_name
%==========================================================================

%==========================================================================
% STL binary file format
%=======================
% Binary STL files have an 84 byte header followed by 50-byte records, each
% describing a single facet of the mesh.  Technically each facet could be
% any 2D shape, but that would screw up the 50-byte-per-facet structure, so
% in practice only triangular facets are used.  The present code ONLY works
% for meshes composed of triangular facets.
%
% HEADER:
% 80 bytes:  Header text
% 4 bytes:   (int) The number of facets in the STL mesh
%
% DATA:
% 4 bytes:  (float) normal x
% 4 bytes:  (float) normal y
% 4 bytes:  (float) normal z
% 4 bytes:  (float) vertex1 x
% 4 bytes:  (float) vertex1 y
% 4 bytes:  (float) vertex1 z
% 4 bytes:  (float) vertex2 x
% 4 bytes:  (float) vertex2 y
% 4 bytes:  (float) vertex2 z
% 4 bytes:  (float) vertex3 x
% 4 bytes:  (float) vertex3 y
% 4 bytes:  (float) vertex3 z
% 2 bytes:  Padding to make the data for each facet 50-bytes in length
%   ...and repeat for next facet... 
%==========================================================================

if nargin==2
  stlFORMAT = lower(varargin{1});
else
  stlFORMAT = 'auto';
end

%If necessary, identify whether the STL is ascii or binary:
if strcmp(stlFORMAT,'ascii')==0 && strcmp(stlFORMAT,'binary')==0
  stlFORMAT = IDENTIFY_stl_format(stlFILENAME);
end

%Load the STL file:
if strcmp(stlFORMAT,'ascii')
  [coordVERTICES,coordNORMALS,stlNAME] = READ_stlascii(stlFILENAME);
elseif strcmp(stlFORMAT,'binary')
  [coordVERTICES,coordNORMALS] = READ_stlbinary(stlFILENAME);
  stlNAME = 'unnamed_object';
end %if

%Prepare the output arguments
if nargout == 2
  varargout(1) = {coordNORMALS};
elseif nargout == 3
  varargout(1) = {coordNORMALS};
  varargout(2) = {stlNAME};
end

end %function



%==========================================================================
function [stlFORMAT] = IDENTIFY_stl_format(stlFILENAME)
% IDENTIFY_stl_format  Test whether an stl file is ascii or binary

% Open the file:
fidIN = fopen(stlFILENAME);

% Check the file size first, since binary files MUST have a size of 84+(50*n)
fseek(fidIN,0,1);         % Go to the end of the file
fidSIZE = ftell(fidIN);   % Check the size of the file

if rem(fidSIZE-84,50) > 0
    
  stlFORMAT = 'ascii';

else

  % Files with a size of 84+(50*n), might be either ascii or binary...
    
  % Read first 80 characters of the file.
  % For an ASCII file, the data should begin immediately (give or take a few
  % blank lines or spaces) and the first word must be 'solid'.
  % For a binary file, the first 80 characters contains the header.
  % It is bad practice to begin the header of a binary file with the word
  % 'solid', so it can be used to identify whether the file is ASCII or
  % binary.
  fseek(fidIN,0,-1);        % Go to the start of the file
  firsteighty = char(fread(fidIN,80,'uchar')');

  % Trim leading and trailing spaces:
  firsteighty = strtrim(firsteighty);

  % Take the first five remaining characters, and check if these are 'solid':
  firstfive = firsteighty(1:min(5,length(firsteighty)));

  % Double check by reading the last 80 characters of the file.
  % For an ASCII file, the data should end (give or take a few
  % blank lines or spaces) with 'endsolid <object_name>'.
  % If the last 80 characters contains the word 'endsolid' then this
  % confirms that the file is indeed ASCII.
  if strcmp(firstfive,'solid')
  
    fseek(fidIN,-80,1);     % Go to the end of the file minus 80 characters
    lasteighty = char(fread(fidIN,80,'uchar')');
  
    if findstr(lasteighty,'endsolid')
      stlFORMAT = 'ascii';
    else
      stlFORMAT = 'binary';
    end
  
  else
    stlFORMAT = 'binary';
  end
  
end

% Close the file
fclose(fidIN);

end %function
%==========================================================================



%==========================================================================
function [coordVERTICES,coordNORMALS,stlNAME] = READ_stlascii(stlFILENAME)
% READ_stlascii  Read mesh data in the form of an ascii <*.stl> file

% Read the ascii STL file
fidIN = fopen(stlFILENAME);
fidCONTENTcell = textscan(fidIN,'%s','delimiter','\n');                  %Read all the file
fidCONTENT = fidCONTENTcell{:}(logical(~strcmp(fidCONTENTcell{:},'')));  %Remove all blank lines
fclose(fidIN);

% Read the STL name
if nargout == 3
  line1 = char(fidCONTENT(1));
  if (size(line1,2) >= 7)
    stlNAME = line1(7:end);
  else
    stlNAME = 'unnamed_object'; 
  end
end

% Read the vector normals
if nargout >= 2
  stringNORMALS = char(fidCONTENT(logical(strncmp(fidCONTENT,'facet normal',12))));
  coordNORMALS  = str2num(stringNORMALS(:,13:end));
end

% Read the vertex coordinates
facetTOTAL       = sum(strcmp(fidCONTENT,'endfacet'));
stringVERTICES   = char(fidCONTENT(logical(strncmp(fidCONTENT,'vertex',6))));
coordVERTICESall = str2num(stringVERTICES(:,7:end));
cotemp           = zeros(3,facetTOTAL,3);
cotemp(:)        = coordVERTICESall;
coordVERTICES    = shiftdim(cotemp,1);

end %function
%==========================================================================



%==========================================================================
function [coordVERTICES,coordNORMALS] = READ_stlbinary(stlFILENAME)
% READ_stlbinary  Read mesh data in the form of an binary <*.stl> file

% Open the binary STL file
fidIN = fopen(stlFILENAME);

% Read the header
fseek(fidIN,80,-1);                   % Move to the last 4 bytes of the header
facetcount = fread(fidIN,1,'int32');  % Read the number of facets in the stl file

% Initialise arrays into which the STL data will be loaded:
coordNORMALS  = zeros(facetcount,3);
coordVERTICES = zeros(facetcount,3,3);

% Read the data for each facet:
for loopF = 1:facetcount
  
  tempIN = fread(fidIN,3*4,'float');
  
  coordNORMALS(loopF,1:3)    = tempIN(1:3);    % x,y,z components of the facet's normal vector
  coordVERTICES(loopF,1:3,1) = tempIN(4:6);    % x,y,z coordinates of vertex 1
  coordVERTICES(loopF,1:3,2) = tempIN(7:9);    % x,y,z coordinates of vertex 2
  coordVERTICES(loopF,1:3,3) = tempIN(10:12);  % x,y,z coordinates of vertex 3
  
  fread(fidIN,1,'int16');   % Move to the start of the next facet.  Using fread is much quicker than using fseek(fidIN,2,0); 

end %for

% Close the binary STL file
fclose(fidIN);

end %function
%==========================================================================
