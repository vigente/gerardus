%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spherical Harmonic Modeling and Analysis Toolkit (SPHARM-MAT) is a 3D 
% shape modeling and analysis toolkit. 
% It is a software package developed at Shenlab in Center for Neuroimaging, 
% Indiana University (SpharmMat@gmail.com, http://www.iupui.edu/~shenlab/)
% It is available to the scientific community as copyright freeware 
% under the terms of the GNU General Public Licence.
% 
% Copyright 2009, 2010, ShenLab, Center for Neuroimaging, Indiana University
% 
% This file is part of SPHARM-MAT.
% 
% SPHARM-MAT is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% SPHARM-MAT is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with SPHARM-MAT. If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fout, vout, cout] = stlreadBIN(filename)
% This function reads an STL file in binary format into matrixes X, Y and
% Z, and C.  C is optional and contains color rgb data in 5 bits.  
%
% USAGE: [x, y, z, c] = stlread(filename);
%
% To plot use patch(x,y,z,c), or patch(x,y,z)
%
% Written by Doron Harlev

if nargout>4
    error('Too many output arguments')
end
use_color=(nargout==4);

fid=fopen(filename, 'r'); %Open the file, assumes STL Binary format.
if fid == -1 
    error('File could not be opened, check name or path.')
end

ftitle=fread(fid,80,'uchar=>schar'); % Read file title
num_facet=fread(fid,1,'int32'); % Read number of Facets

%fprintf('\nTitle: %s\n', char(ftitle'));
%fprintf('Num Facets: %d\n', num_facet);

% Preallocate memory to save running time
x=zeros(3,num_facet); y=zeros(3,num_facet); z=zeros(3,num_facet);
if use_color
    c=uint8(zeros(3,num_facet));
end

[pa,na,ex,ve] = fileparts(filename);

str=sprintf('Reading %s ...', na);

h = waitbar(0,str);
for i=1:num_facet,
    norm=fread(fid,3,'float32'); % normal coordinates, ignored for now
    ver1=fread(fid,3,'float32'); % vertex 1
    ver2=fread(fid,3,'float32'); % vertex 2
    ver3=fread(fid,3,'float32'); % vertex 3
    col=fread(fid,1,'uint16'); % color bytes
    if (bitget(col,16)==1 & use_color)
        r=bitshift(bitand(2^16-1, col),-10);
        g=bitshift(bitand(2^11-1, col),-5);
        b=bitand(2^6-1, col);
        c(:,i)=[r; g; b];
    end
    x(:,i)=[ver1(1); ver2(1); ver3(1)]; % convert to matlab "patch" compatible format
    y(:,i)=[ver1(2); ver2(2); ver3(2)];
    z(:,i)=[ver1(3); ver2(3); ver3(3)];
    if mod(i,floor(num_facet/10))==0
        waitbar(i/num_facet,h);
    end
end

if use_color
    cout(1)={c};
else
    cout = [];
end
fclose(fid);
close(h);

h = figure(1);
p = patch(x,y,z,z);

vout = get(p,'Vertices');
fout = get(p,'Faces');

close(h);

return;
% For more information http://rpdrc.ic.polyu.edu.hk/old_files/stl_binary_format.htm
