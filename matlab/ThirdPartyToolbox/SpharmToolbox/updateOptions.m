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

function updateOptions(handles, confs)

for k = 1:length(confs.vars)
    eval(sprintf('vHandle=confs.%s;',confs.vars{k}));
    strV = sprintf('handles.v%d', k);
    if ~isempty(vHandle)
        if confs.args(k) < 10
            set(eval(strV), 'String', num2str(vHandle));
        else
            vHandle = strtrim(vHandle);

            sty = get(eval(strV),'Style');
            if strcmpi(sty,'popupmenu')
                contentStr = get(eval(strV),'String');
                idx = find(strcmpi(vHandle, strtrim(contentStr))==1);
                set(eval(strV), 'Value', idx);
            else
                set(eval(strV), 'String', vHandle);
            end
        end
    end
end

if ~isempty(handles.userdata.inObjs)
    set(handles.IOL_tag, 'String', handles.userdata.inObjs);
end
if ~isempty(handles.userdata.inObjs2)
    set(handles.IOL2_tag, 'String', handles.userdata.inObjs2);
end

return;