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

% ============================================
% gipl_read_volume.m
% function for reading volume of Guys Image Processing Lab (Gipl) volume
% file
% INPUT: info, which is head information readed from gipl file. 
% 
% OUTPUT: V, which store the datum of the gipl file. 
% 
% Wan Jing
% 11/06/2008 - create
%============================================
function [V] = gipl_read_volume(info)

if(~isstruct(info)) info=gipl_read_header(info); end

% Open gipl file
f=fopen(getfield(info,'filename'),'rb','ieee-be');

  % Seek volume data start
  if(info.image_type==1), voxelbits=1; end
  if(info.image_type==7||info.image_type==8), voxelbits=8; end
  if(info.image_type==15||info.image_type==16), voxelbits=16; end
  if(info.image_type==31||info.image_type==32||info.image_type==64), voxelbits=32; end
  if(info.image_type==65), voxelbits=64; end
 

  % Read Volume data
  tmp=getfield(info,'sizes');
  volsize(1:3)=tmp(1:3);
  
  datasize=prod(tmp(1:3))*(voxelbits/8);
  fsize=getfield(info,'filesize');
  fseek(f,fsize-datasize,'bof');
  
  if(info.image_type==1), V = logical(fread(f,datasize,'bit1')); end
  if(info.image_type==7), V = int8(fread(f,datasize,'char')); end
  if(info.image_type==8), V = uint8(fread(f,datasize,'uchar')); end
  if(info.image_type==15), V = int16(fread(f,datasize,'short')); end 
  if(info.image_type==16), V = uint16(fread(f,datasize,'ushort')); end
  if(info.image_type==31), V = uint32(fread(f,datasize,'uint')); end
  if(info.image_type==32), V = int32(fread(f,datasize,'int')); end
  if(info.image_type==64), V = single(fread(f,datasize,'float')); end 
  if(info.image_type==65), V = double(fread(f,datasize,'double')); end 

fclose(f);

% Reshape the volume data to the right dimensions
V = reshape(V,volsize);

return;