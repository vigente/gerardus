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

function [inObjs1, inObjs2] = LoadObjList(inputFile)

fin = fopen(inputFile, 'rt');
inObjs1={};
inObjs2={};

[numOut, vals] = Get_Vals(fin, 'inputs', 10^6);
if numOut >0
    inObjs1 = vals;
%    set(handles.IOL_tag, 'String', handles.userdata.inObjs);
end

[numOut, vals] = Get_Vals(fin, 'inputs2', 10^6);
if numOut >0
    inObjs2 = vals;
%    set(handles.IOL_tag, 'String', handles.userdata.inObjs);
end

fclose(fin);
    
return;


function [numOut, vals] = Get_Vals(fInput, Keyword, NumVals)

%global userdata;

numOut = 0;
vals = {};

fseek(fInput, 0, 'bof');
delim = ' ,';

while 1
    tline = fgetl(fInput);
    if ~ischar(tline),   break,   end;
    [key, remain] = strtok(tline, delim);
    if strcmpi(strtrim(key), strtrim(Keyword))
        while numOut < NumVals
            tline = fgetl(fInput);
            if ~ischar(tline),   break,   end;
            numOut = numOut + 1;
            vals(numOut) = {deblank(tline)};
        end
        break;
    end
end    

return;