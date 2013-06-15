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

function handles = dispOptions(handles)

handles.userdata.IOT_tag_pos = get(handles.IOT_tag, 'Position');
handles.userdata.IO_tag_pos = get(handles.IO_tag, 'Position');
handles.userdata.IOL_tag_pos = get(handles.IOL_tag, 'Position');

pos = get(handles.v1,'Position');
height = (pos(4)*1.2);

switch handles.userdata.keyword.selected
    case 'ParamCALD'
        for i = 1:length(handles.userdata.ParamCALD.vars)
            strVT = sprintf('handles.vt%d',i);
            strV = sprintf('handles.v%d', i);
            strVC = sprintf('handles.vc%d', i);
            set(eval(strVT), 'Visible', 'on', 'String', handles.userdata.ParamCALD.vars{i});

            if strcmp(handles.userdata.ParamCALD.vars{i},'Mesh')
                set(eval(strV), 'Visible', 'on', 'String', handles.userdata.mesh,'Style','popupmenu','Value',1);
                contentStr = get(eval(strV),'String')
                eval(sprintf('handles.userdata.ParamCALD.%s=contentStr{valsMenu};',handles.userdata.ParamCALD.vars{i}));
            elseif strcmp(handles.userdata.ParamCALD.vars{i},'t_major')
                set(eval(strV), 'Visible', 'on', 'String', {'x';'y'},'Style','popupmenu');  
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.ParamCALD.%s=contentStr{valsMenu};',handles.userdata.ParamCALD.vars{i}));
            elseif strcmp(handles.userdata.ParamCALD.vars{i},'SelectDiagonal');
                set(eval(strV), 'Visible', 'on', 'String', {'ShortDiag';'LongDiag'},'Style','popupmenu');
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.ParamCALD.%s=contentStr{valsMenu};',handles.userdata.ParamCALD.vars{i}));
            else
                cValStr = ['handles.userdata.ParamCALD.' handles.userdata.ParamCALD.vars{i}];
                cVal = eval(cValStr);
                if handles.userdata.ParamCALD.args(i) < 10
                    set(eval(strV), 'Visible', 'on', 'String', num2str(cVal));
                elseif handles.userdata.ParamCALD.args(i) >= 10
                    set(eval(strV), 'Visible', 'on', 'String', cVal);
                end                
%                 set(eval(strV), 'Visible', 'on', 'String', handles.userdata.ParamCALD.default{i});
            end

            if (handles.userdata.ParamCALD.args(i) >= 100)
                set(eval(strVC), 'Visible', 'on');
            else
                set(eval(strVC), 'Visible', 'off');
            end
        end
        for j = length(handles.userdata.ParamCALD.vars)+1:10
            strVT = sprintf('handles.vt%d',j);
            strV = sprintf('handles.v%d', j);
            strVC = sprintf('handles.vc%d', j);
            set(eval(strVT), 'Visible', 'off', 'String', '');
            set(eval(strV), 'Visible', 'off', 'String', '');
            set(eval(strVC), 'Visible', 'off');
        end
        set(handles.IOT_tag, 'Position', handles.userdata.IOT_tag_pos+[0 (10-length(handles.userdata.ParamCALD.vars))*height 0 0]);
        set(handles.IO_tag, 'Position', handles.userdata.IO_tag_pos+[0 (10-length(handles.userdata.ParamCALD.vars))*height 0 0]);
        set(handles.IOL_tag, 'Position', handles.userdata.IOL_tag_pos+[0 0 0 (10-length(handles.userdata.ParamCALD.vars))*height]);        
        implemented = 1;
    case 'ParamPDM'
        for i = 1:length(handles.userdata.ParamPDM.vars)
            strVT = sprintf('handles.vt%d',i);
            strV = sprintf('handles.v%d', i);
            strVC = sprintf('handles.vc%d', i);
            set(eval(strVT), 'Visible', 'on', 'String', handles.userdata.ParamPDM.vars{i});

            if strcmp(handles.userdata.ParamPDM.vars{i},'Mesh');
                set(eval(strV), 'Visible', 'on', 'String', handles.userdata.mesh,'Style','popupmenu');
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.ParamPDM.%s=contentStr{valsMenu};',handles.userdata.ParamPDM.vars{i}));
            else
                cValStr = ['handles.userdata.ParamPDM.' handles.userdata.ParamPDM.vars{i}];
                cVal = eval(cValStr);
                if handles.userdata.ParamPDM.args(i) < 10
                    set(eval(strV), 'Visible', 'on', 'String', num2str(cVal));
                elseif handles.userdata.ParamPDM.args(i) >= 10
                    set(eval(strV), 'Visible', 'on', 'String', cVal);
                end
%                 set(eval(strV), 'Visible', 'on', 'String', handles.userdata.ParamPDM.default{i});
            end

            if (handles.userdata.ParamPDM.args(i) >= 100)
                set(eval(strVC), 'Visible', 'on');
            else
                set(eval(strVC), 'Visible', 'off');
            end
        end
        for j = length(handles.userdata.ParamPDM.vars)+1:10
            strVT = sprintf('handles.vt%d',j);
            strV = sprintf('handles.v%d', j);
            set(eval(strVT), 'Visible', 'off', 'String', '');
            set(eval(strV), 'Visible', 'off', 'String', '');
            strVC = sprintf('handles.vc%d', j);
            set(eval(strVC), 'Visible', 'off');
        end
        set(handles.IOT_tag, 'Position', handles.userdata.IOT_tag_pos+[0 (10-length(handles.userdata.ParamPDM.vars))*height 0 0]);
        set(handles.IO_tag, 'Position', handles.userdata.IO_tag_pos+[0 (10-length(handles.userdata.ParamPDM.vars))*height 0 0]);
        set(handles.IOL_tag, 'Position', handles.userdata.IOL_tag_pos+[0 0 0 (10-length(handles.userdata.ParamPDM.vars))*height]);        
        implemented = 1;
    case 'ExpLSF'
        for i = 1:length(handles.userdata.ExpLSF.vars)
            strVT = sprintf('handles.vt%d',i);
            strV = sprintf('handles.v%d', i);
            strVC = sprintf('handles.vc%d', i);
            set(eval(strVT), 'Visible', 'on', 'String', handles.userdata.ExpLSF.vars{i});

            if strcmp(handles.userdata.ExpLSF.vars{i},'Mesh');
                set(eval(strV), 'Visible', 'on', 'String', handles.userdata.mesh,'Style','popupmenu');
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.ExpLSF.%s=contentStr{valsMenu};',handles.userdata.ExpLSF.vars{i}));
            else
                cValStr = ['handles.userdata.ExpLSF.' handles.userdata.ExpLSF.vars{i}];
                cVal = eval(cValStr);
                if handles.userdata.ExpLSF.args(i) < 10
                    set(eval(strV), 'Visible', 'on', 'String', num2str(cVal));
                elseif handles.userdata.ExpLSF.args(i) >= 10
                    set(eval(strV), 'Visible', 'on', 'String', cVal);
                end
%                 set(eval(strV), 'Visible', 'on', 'String', handles.userdata.ExpLSF.default{i});
            end

            if (handles.userdata.ExpLSF.args(i) >= 100)
                set(eval(strVC), 'Visible', 'on');
            else
                set(eval(strVC), 'Visible', 'off');
            end
        end
        for j = length(handles.userdata.ExpLSF.vars)+1:10
            strVT = sprintf('handles.vt%d',j);
            strV = sprintf('handles.v%d', j);
            set(eval(strVT), 'Visible', 'off', 'String', '');
            set(eval(strV), 'Visible', 'off', 'String', '');
            strVC = sprintf('handles.vc%d', j);
            set(eval(strVC), 'Visible', 'off');
        end
        set(handles.IOT_tag, 'Position', handles.userdata.IOT_tag_pos+[0 (10-length(handles.userdata.ExpLSF.vars))*height 0 0]);
        set(handles.IO_tag, 'Position', handles.userdata.IO_tag_pos+[0 (10-length(handles.userdata.ExpLSF.vars))*height 0 0]);
        set(handles.IOL_tag, 'Position', handles.userdata.IOL_tag_pos+[0 0 0 (10-length(handles.userdata.ExpLSF.vars))*height]);
        implemented = 1;
    case 'AligSHREC'
        for i = 1:length(handles.userdata.AligSHREC.vars)
            strVT = sprintf('handles.vt%d',i);
            strV = sprintf('handles.v%d', i);
            strVC = sprintf('handles.vc%d', i);
            set(eval(strVT), 'Visible', 'on', 'String', handles.userdata.AligSHREC.vars{i});
            if strcmp(handles.userdata.AligSHREC.vars{i},'Mesh')
                set(eval(strV), 'Visible', 'on', 'String', handles.userdata.mesh,'Style','popupmenu');
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.AligSHREC.%s=contentStr{valsMenu};',handles.userdata.AligSHREC.vars{i}));
            elseif strcmp(handles.userdata.AligSHREC.vars{i},'NormalizeSize')
                set(eval(strV), 'Visible', 'on', 'String', {'No';'Yes'},'Style','popupmenu');   
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.AligSHREC.%s=contentStr{valsMenu};',handles.userdata.AligSHREC.vars{i}));
            else
                cValStr = ['handles.userdata.AligSHREC.' handles.userdata.AligSHREC.vars{i}];
                cVal = eval(cValStr);
                if handles.userdata.AligSHREC.args(i) < 10
                    set(eval(strV), 'Visible', 'on', 'String', num2str(cVal));
                elseif handles.userdata.AligSHREC.args(i) >= 10
                    set(eval(strV), 'Visible', 'on', 'String', cVal);
                end
%                 set(eval(strV), 'Visible', 'on', 'String', handles.userdata.AligSHREC.default{i});
            end

            if (handles.userdata.AligSHREC.args(i) >= 100)
                set(eval(strVC), 'Visible', 'on');
            else
                set(eval(strVC), 'Visible', 'off');
            end
        end
%        set(handles.vc1, 'Visible', 'on');
        for j = length(handles.userdata.AligSHREC.vars)+1:10
            strVT = sprintf('handles.vt%d',j);
            strV = sprintf('handles.v%d', j);
            set(eval(strVT), 'Visible', 'off', 'String', '');
            set(eval(strV), 'Visible', 'off', 'String', '');
            strVC = sprintf('handles.vc%d', j);
            set(eval(strVC), 'Visible', 'off');
        end
        set(handles.IOT_tag, 'Position', handles.userdata.IOT_tag_pos+[0 (10-length(handles.userdata.AligSHREC.vars))*height 0 0]);
        set(handles.IO_tag, 'Position', handles.userdata.IO_tag_pos+[0 (10-length(handles.userdata.AligSHREC.vars))*height 0 0]);
        set(handles.IOL_tag, 'Position', handles.userdata.IOL_tag_pos+[0 0 0 (10-length(handles.userdata.AligSHREC.vars))*height]);
        implemented = 1;
    case 'AligFOE'
        for i = 1:length(handles.userdata.AligFOE.vars)
            strVT = sprintf('handles.vt%d',i);
            strV = sprintf('handles.v%d', i);
            strVC = sprintf('handles.vc%d', i);
            set(eval(strVT), 'Visible', 'on', 'String', handles.userdata.AligFOE.vars{i});
            if strcmp(handles.userdata.AligFOE.vars{i},'Mesh')
                set(eval(strV), 'Visible', 'on', 'String', handles.userdata.mesh,'Style','popupmenu');
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.AligFOE.%s=contentStr{valsMenu};',handles.userdata.AligFOE.vars{i}));
            elseif strcmp(handles.userdata.AligFOE.vars{i},'CPoint')
                set(eval(strV), 'Visible', 'on', 'String', {'x';'y';'z'},'Style','popupmenu','Value',2);
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.AligFOE.%s=contentStr{valsMenu};',handles.userdata.AligFOE.vars{i}));
            elseif strcmp(handles.userdata.AligFOE.vars{i},'NPole')
                set(eval(strV), 'Visible', 'on', 'String', {'x';'y';'z'},'Style','popupmenu','Value',3);
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.AligFOE.%s=contentStr{valsMenu};',handles.userdata.AligFOE.vars{i}));
            else
                cValStr = ['handles.userdata.AligFOE.' handles.userdata.AligFOE.vars{i}];
                cVal = eval(cValStr);
                if handles.userdata.AligFOE.args(i) < 10
                    set(eval(strV), 'Visible', 'on', 'String', num2str(cVal));
                elseif handles.userdata.AligFOE.args(i) >= 10
                    set(eval(strV), 'Visible', 'on', 'String', cVal);
                end
%                 set(eval(strV), 'Visible', 'on', 'String', handles.userdata.AligFOE.default{i});
            end

            if (handles.userdata.AligFOE.args(i) >= 100)
                set(eval(strVC), 'Visible', 'on');
            else
                set(eval(strVC), 'Visible', 'off');
            end
        end
        for j = length(handles.userdata.AligFOE.vars)+1:10
            strVT = sprintf('handles.vt%d',j);
            strV = sprintf('handles.v%d', j);
            set(eval(strVT), 'Visible', 'off', 'String', '');
            set(eval(strV), 'Visible', 'off', 'String', '');
            strVC = sprintf('handles.vc%d', j);
            set(eval(strVC), 'Visible', 'off');
        end
        set(handles.IOT_tag, 'Position', handles.userdata.IOT_tag_pos+[0 (10-length(handles.userdata.AligFOE.vars))*height 0 0]);
        set(handles.IO_tag, 'Position', handles.userdata.IO_tag_pos+[0 (10-length(handles.userdata.AligFOE.vars))*height 0 0]);
        set(handles.IOL_tag, 'Position', handles.userdata.IOL_tag_pos+[0 0 0 (10-length(handles.userdata.AligFOE.vars))*height]);
        implemented = 1;
    case 'ExpAlig'
        for i = 1:length(handles.userdata.ExpAlig.vars)
            strVT = sprintf('handles.vt%d',i);
            strV = sprintf('handles.v%d', i);
            strVC = sprintf('handles.vc%d', i);
            set(eval(strVT), 'Visible', 'on', 'String', handles.userdata.ExpAlig.vars{i});
            if strcmp(handles.userdata.ExpAlig.vars{i},'Mesh');
                set(eval(strV), 'Visible', 'on', 'String', handles.userdata.mesh,'Style','popupmenu');
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.ExpAlig.%s=contentStr{valsMenu};',handles.userdata.ExpAlig.vars{i}));
            else
                cValStr = ['handles.userdata.ExpAlig.' handles.userdata.ExpAlig.vars{i}];
                cVal = eval(cValStr);
                if handles.userdata.ExpAlig.args(i) < 10
                    set(eval(strV), 'Visible', 'on', 'String', num2str(cVal));
                elseif handles.userdata.ExpAlig.args(i) >= 10
                    set(eval(strV), 'Visible', 'on', 'String', cVal);
                end
%                 set(eval(strV), 'Visible', 'on', 'String', handles.userdata.ExpAlig.default{i});
            end

            if (handles.userdata.ExpAlig.args(i) >= 100)
                set(eval(strVC), 'Visible', 'on');
            else
                set(eval(strVC), 'Visible', 'off');
            end
        end
        for j = length(handles.userdata.ExpAlig.vars)+1:10
            strVT = sprintf('handles.vt%d',j);
            strV = sprintf('handles.v%d', j);
            set(eval(strVT), 'Visible', 'off', 'String', '');
            set(eval(strV), 'Visible', 'off', 'String', '');
            strVC = sprintf('handles.vc%d', j);
            set(eval(strVC), 'Visible', 'off');
        end
        set(handles.IOT_tag, 'Position', handles.userdata.IOT_tag_pos+[0 (10-length(handles.userdata.ExpAlig.vars))*height 0 0]);
        set(handles.IO_tag, 'Position', handles.userdata.IO_tag_pos+[0 (10-length(handles.userdata.ExpAlig.vars))*height 0 0]);
        set(handles.IOL_tag, 'Position', handles.userdata.IOL_tag_pos+[0 0 0 (10-length(handles.userdata.ExpAlig.vars))*height]);
        implemented = 1;
    case 't_map'
        for i = 1:length(handles.userdata.t_map.vars)
            strVT = sprintf('handles.vt%d',i);
            strV = sprintf('handles.v%d', i);
            strVC = sprintf('handles.vc%d', i);
            set(eval(strVT), 'Visible', 'on', 'String', handles.userdata.t_map.vars{i});
            
            if strcmp(handles.userdata.t_map.vars{i},'SampleMesh')
                set(eval(strV), 'Visible', 'on', 'String', handles.userdata.mesh(7:end),'Style','popupmenu');
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.t_map.%s=contentStr{valsMenu};',handles.userdata.t_map.vars{i}));
            elseif strcmp(handles.userdata.t_map.vars{i},'EqualVariance')
                set(eval(strV), 'Visible', 'on', 'String', {'Yes';'No'},'Style','popupmenu');                
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.t_map.%s=contentStr{valsMenu};',handles.userdata.t_map.vars{i}));
            elseif strcmp(handles.userdata.t_map.vars{i},'Signal')
                set(eval(strV), 'Visible', 'on', 'String', {'vl_defm_org';'vl_defm_nrm';'vl_defm_pca';'vl_defm_fld'},'Style','popupmenu');
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.t_map.%s=contentStr{valsMenu};',handles.userdata.t_map.vars{i}));
            else
                cValStr = ['handles.userdata.t_map.' handles.userdata.t_map.vars{i}];
                cVal = eval(cValStr);
                if handles.userdata.t_map.args(i) < 10
                    set(eval(strV), 'Visible', 'on', 'String', num2str(cVal));
                elseif handles.userdata.t_map.args(i) >= 10
                    set(eval(strV), 'Visible', 'on', 'String', cVal);
                end
%                 set(eval(strV), 'Visible', 'on', 'String', handles.userdata.t_map.default{i});
            end

            if (handles.userdata.t_map.args(i) >= 100)
                set(eval(strVC), 'Visible', 'on');
            else
                set(eval(strVC), 'Visible', 'off');
            end
        end
        
        for j = length(handles.userdata.t_map.vars)+1:10
            strVT = sprintf('handles.vt%d',j);
            strV = sprintf('handles.v%d', j);
            set(eval(strVT), 'Visible', 'off', 'String', '');
            set(eval(strV), 'Visible', 'off', 'String', '');
            strVC = sprintf('handles.vc%d', j);
            set(eval(strVC), 'Visible', 'off');
        end
        set(handles.IOT_tag, 'String','Group1','Position', handles.userdata.IOT_tag_pos+[0 (10-length(handles.userdata.t_map.vars))*(height*1.2) 0 0]);
        set(handles.IO_tag, 'Position', handles.userdata.IO_tag_pos+[0 (10-length(handles.userdata.t_map.vars))*(height*1.2) 0 0]);
        set(handles.IOL_tag, 'Position', handles.userdata.IOL_tag_pos+[0 ((10-length(handles.userdata.t_map.vars)))*(height*1.2) 0 0]);

        set(handles.IOT2_tag, 'Visible','on', 'String','Group2','Position', handles.userdata.IOT_tag_pos+[0 -(height*0.1) 0 0]);
        set(handles.IO2_tag, 'Visible','on', 'Position', handles.userdata.IO_tag_pos+[0 -(height*0.1) 0 0]);
        set(handles.IOL2_tag, 'Visible','on', 'Position', handles.userdata.IOL_tag_pos+[0 -(height*0.1) 0 0]);
        implemented = 1;
    case 'PCA'
        for i = 1:length(handles.userdata.PCA.vars)
            strVT = sprintf('handles.vt%d',i);
            strV = sprintf('handles.v%d', i);
            strVC = sprintf('handles.vc%d', i);
            set(eval(strVT), 'Visible', 'on', 'String', handles.userdata.PCA.vars{i});
            
            cValStr = ['handles.userdata.PCA.' handles.userdata.PCA.vars{i}];
            cVal = eval(cValStr);
            if handles.userdata.PCA.args(i) < 10
                set(eval(strV), 'Visible', 'on', 'String', num2str(cVal));
            elseif handles.userdata.PCA.args(i) >= 10
                set(eval(strV), 'Visible', 'on', 'String', cVal);
            end
            
%             set(eval(strV), 'Visible', 'on', 'String', handles.userdata.PCA.default{i});

            if (handles.userdata.PCA.args(i) >= 100)
                set(eval(strVC), 'Visible', 'on');
            else
                set(eval(strVC), 'Visible', 'off');
            end
        end
        for j = length(handles.userdata.PCA.vars)+1:10
            strVT = sprintf('handles.vt%d',j);
            strV = sprintf('handles.v%d', j);
            set(eval(strVT), 'Visible', 'off', 'String', '');
            set(eval(strV), 'Visible', 'off', 'String', '');
            strVC = sprintf('handles.vc%d', j);
            set(eval(strVC), 'Visible', 'off');
        end
        set(handles.IOT_tag, 'Position', handles.userdata.IOT_tag_pos+[0 (10-length(handles.userdata.PCA.vars))*height 0 0]);
        set(handles.IO_tag, 'Position', handles.userdata.IO_tag_pos+[0 (10-length(handles.userdata.PCA.vars))*height 0 0]);
        set(handles.IOL_tag, 'Position', handles.userdata.IOL_tag_pos+[0 0 0 (10-length(handles.userdata.PCA.vars))*height]);
        implemented = 1;
    case 'res_t_map'
        for i = 1:length(handles.userdata.res_t_map.vars)
            strVT = sprintf('handles.vt%d',i);
            strV = sprintf('handles.v%d', i);
            strVC = sprintf('handles.vc%d', i);
            set(eval(strVT), 'Visible', 'on', 'String', handles.userdata.res_t_map.vars{i});
            if strcmp(handles.userdata.res_t_map.vars{i},'Colormap')
                set(eval(strV), 'Visible', 'on', 'String', handles.userdata.colormap,'Style','popupmenu');
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.res_t_map.%s=contentStr{valsMenu};',handles.userdata.res_t_map.vars{i}));
            elseif strcmp(handles.userdata.res_t_map.vars{i},'Overlay')
                set(eval(strV), 'Visible', 'on', 'String', {'p-value';'t-map'},'Style','popupmenu');  
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.res_t_map.%s=contentStr{valsMenu};',handles.userdata.res_t_map.vars{i}));
            else
                cValStr = ['handles.userdata.res_t_map.' handles.userdata.res_t_map.vars{i}];
                cVal = eval(cValStr);
                if handles.userdata.res_t_map.args(i) < 10
                    set(eval(strV), 'Visible', 'on', 'String', num2str(cVal));
                elseif handles.userdata.res_t_map.args(i) >= 10
                    set(eval(strV), 'Visible', 'on', 'String', cVal);
                end
%                 set(eval(strV), 'Visible', 'on', 'String', handles.userdata.res_t_map.default{i});
            end

            if (handles.userdata.res_t_map.args(i) >= 100)
                set(eval(strVC), 'Visible', 'on');
            else
                set(eval(strVC), 'Visible', 'off');
            end
        end
        for j = length(handles.userdata.res_t_map.vars)+1:10
            strVT = sprintf('handles.vt%d',j);
            strV = sprintf('handles.v%d', j);
            set(eval(strVT), 'Visible', 'off', 'String', '');
            set(eval(strV), 'Visible', 'off', 'String', '');
            strVC = sprintf('handles.vc%d', j);
            set(eval(strVC), 'Visible', 'off');
        end
        set(handles.IOT_tag, 'Position', handles.userdata.IOT_tag_pos+[0 (10-length(handles.userdata.res_t_map.vars))*height 0 0]);
        set(handles.IO_tag, 'Position', handles.userdata.IO_tag_pos+[0 (10-length(handles.userdata.res_t_map.vars))*height 0 0]);
        set(handles.IOL_tag, 'Position', handles.userdata.IOL_tag_pos+[0 0 0 (10-length(handles.userdata.res_t_map.vars))*height]);
        implemented = 1;
    case 'res_PCA'
        for i = 1:length(handles.userdata.res_PCA.vars)
            strVT = sprintf('handles.vt%d',i);
            strV = sprintf('handles.v%d', i);
            strVC = sprintf('handles.vc%d', i);
            set(eval(strVT), 'Visible', 'on', 'String', handles.userdata.res_PCA.vars{i});
            if strcmp(handles.userdata.res_PCA.vars{i},'Mesh');
                set(eval(strV), 'Visible', 'on', 'String', handles.userdata.mesh(2:end),'Style','popupmenu');
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.res_PCA.%s=contentStr{valsMenu};',handles.userdata.res_PCA.vars{i}));
            else
                cValStr = ['handles.userdata.res_PCA.' handles.userdata.res_PCA.vars{i}];
                cVal = eval(cValStr);
                if handles.userdata.res_PCA.args(i) < 10
                    set(eval(strV), 'Visible', 'on', 'String', num2str(cVal));
                elseif handles.userdata.res_PCA.args(i) >= 10
                    set(eval(strV), 'Visible', 'on', 'String', cVal);
                end
%                 set(eval(strV), 'Visible', 'on', 'String', handles.userdata.res_PCA.default{i});
            end

            if (handles.userdata.res_PCA.args(i) >= 100)
                set(eval(strVC), 'Visible', 'on');
            else
                set(eval(strVC), 'Visible', 'off');
            end
        end
        for j = length(handles.userdata.res_PCA.vars)+1:10
            strVT = sprintf('handles.vt%d',j);
            strV = sprintf('handles.v%d', j);
            set(eval(strVT), 'Visible', 'off', 'String', '');
            set(eval(strV), 'Visible', 'off', 'String', '');
            strVC = sprintf('handles.vc%d', j);
            set(eval(strVC), 'Visible', 'off');
        end
        set(handles.IOT_tag, 'Position', handles.userdata.IOT_tag_pos+[0 (10-length(handles.userdata.res_PCA.vars))*height 0 0]);
        set(handles.IO_tag, 'Position', handles.userdata.IO_tag_pos+[0 (10-length(handles.userdata.res_PCA.vars))*height 0 0]);
        set(handles.IOL_tag, 'Position', handles.userdata.IOL_tag_pos+[0 0 0 (10-length(handles.userdata.res_PCA.vars))*height]);
        implemented = 1;
    case 'DisplayObjs'
        for i = 1:length(handles.userdata.DisplayObjs.vars)
            strVT = sprintf('handles.vt%d',i);
            strV = sprintf('handles.v%d', i);
            strVC = sprintf('handles.vc%d', i);
            set(eval(strVT), 'Visible', 'on', 'String', handles.userdata.DisplayObjs.vars{i});
            if strcmp(handles.userdata.DisplayObjs.vars{i},'Mesh')
                set(eval(strV), 'Visible', 'on', 'String', handles.userdata.mesh,'Style','popupmenu');
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.DisplayObjs.%s=contentStr{valsMenu};',handles.userdata.DisplayObjs.vars{i}));
            elseif strcmp(handles.userdata.DisplayObjs.vars{i}, 'Shade')
                set(eval(strV), 'Visible', 'on', 'String', handles.userdata.shade,'Style','popupmenu');                
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.DisplayObjs.%s=contentStr{valsMenu};',handles.userdata.DisplayObjs.vars{i}));
            elseif strcmp(handles.userdata.DisplayObjs.vars{i}, 'Space')
                set(eval(strV), 'Visible', 'on', 'String', handles.userdata.space,'Style','popupmenu');
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.DisplayObjs.%s=contentStr{valsMenu};',handles.userdata.DisplayObjs.vars{i}));
            elseif strcmp(handles.userdata.DisplayObjs.vars{i}, 'Overlay')
                set(eval(strV), 'Visible', 'on', 'String', handles.userdata.overlay,'Style','popupmenu');
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.DisplayObjs.%s=contentStr{valsMenu};',handles.userdata.DisplayObjs.vars{i}));
            elseif strcmp(handles.userdata.DisplayObjs.vars{i}, 'Export')
                set(eval(strV), 'Visible', 'on', 'String', handles.userdata.export,'Style','popupmenu');
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.DisplayObjs.%s=contentStr{valsMenu};',handles.userdata.DisplayObjs.vars{i}));
%             elseif strcmp(handles.userdata.DisplayObjs.vars{i}, 'Degree')
%                 set(eval(strV), 'Visible', 'on', 'String', num2str(handles.userdata.DisplayObjs.Degree));
            else
                cValStr = ['handles.userdata.DisplayObjs.' handles.userdata.DisplayObjs.vars{i}];
                cVal = eval(cValStr);
                if handles.userdata.DisplayObjs.args(i) < 10
                    set(eval(strV), 'Visible', 'on', 'String', num2str(cVal));
                elseif handles.userdata.DisplayObjs.args(i) >= 10
                    set(eval(strV), 'Visible', 'on', 'String', cVal);
                end                
%                 set(eval(strV), 'Visible', 'on', 'String', handles.userdata.DisplayObjs.default{i});
            end

            if (handles.userdata.DisplayObjs.args(i) >= 100)
                set(eval(strVC), 'Visible', 'on');
            else
                set(eval(strVC), 'Visible', 'off');
            end
        end
        for j = length(handles.userdata.DisplayObjs.vars)+1:10
            strVT = sprintf('handles.vt%d',j);
            strV = sprintf('handles.v%d', j);
            set(eval(strVT), 'Visible', 'off', 'String', '');
            set(eval(strV), 'Visible', 'off', 'String', '');
            strVC = sprintf('handles.vc%d', j);
            set(eval(strVC), 'Visible', 'off');
        end
        set(handles.IOT_tag, 'Position', handles.userdata.IOT_tag_pos+[0 (10-length(handles.userdata.DisplayObjs.vars))*height 0 0]);
        set(handles.IO_tag, 'Position', handles.userdata.IO_tag_pos+[0 (10-length(handles.userdata.DisplayObjs.vars))*height 0 0]);
        set(handles.IOL_tag, 'Position', handles.userdata.IOL_tag_pos+[0 0 0 (10-length(handles.userdata.DisplayObjs.vars))*height]);
        implemented = 1;        
    case 'AverageObjs'
        for i = 1:length(handles.userdata.AverageObjs.vars)
            strVT = sprintf('handles.vt%d',i);
            strV = sprintf('handles.v%d', i);
            strVC = sprintf('handles.vc%d', i);
            set(eval(strVT), 'Visible', 'on', 'String', handles.userdata.AverageObjs.vars{i});

            cValStr = ['handles.userdata.AverageObjs.' handles.userdata.AverageObjs.vars{i}];
            cVal = eval(cValStr);
            if handles.userdata.AverageObjs.args(i) < 10
                set(eval(strV), 'Visible', 'on', 'String', num2str(cVal));
            elseif handles.userdata.AverageObjs.args(i) >= 10
                set(eval(strV), 'Visible', 'on', 'String', cVal);
            end

            if (handles.userdata.AverageObjs.args(i) >= 100)
                set(eval(strVC), 'Visible', 'on');
            else
                set(eval(strVC), 'Visible', 'off');
            end
        end
        for j = length(handles.userdata.AverageObjs.vars)+1:10
            strVT = sprintf('handles.vt%d',j);
            strV = sprintf('handles.v%d', j);
            set(eval(strVT), 'Visible', 'off', 'String', '');
            set(eval(strV), 'Visible', 'off', 'String', '');
            strVC = sprintf('handles.vc%d', j);
            set(eval(strVC), 'Visible', 'off');
        end
        set(handles.IOT_tag, 'Position', handles.userdata.IOT_tag_pos+[0 (10-length(handles.userdata.AverageObjs.vars))*height 0 0]);
        set(handles.IO_tag, 'Position', handles.userdata.IO_tag_pos+[0 (10-length(handles.userdata.AverageObjs.vars))*height 0 0]);
        set(handles.IOL_tag, 'Position', handles.userdata.IOL_tag_pos+[0 0 0 (10-length(handles.userdata.AverageObjs.vars))*height]);
        implemented = 1;
    case 'ScaleObjs'
        for i = 1:length(handles.userdata.ScaleObjs.vars)
            strVT = sprintf('handles.vt%d',i);
            strV = sprintf('handles.v%d', i);
            strVC = sprintf('handles.vc%d', i);
            set(eval(strVT), 'Visible', 'on', 'String', handles.userdata.ScaleObjs.vars{i});
            
            
            cValStr = ['handles.userdata.ScaleObjs.' handles.userdata.ScaleObjs.vars{i}];
            cVal = eval(cValStr);
            if handles.userdata.ScaleObjs.args(i) < 10
                set(eval(strV), 'Visible', 'on', 'String', num2str(cVal));
            elseif handles.userdata.ScaleObjs.args(i) >= 10
                set(eval(strV), 'Visible', 'on', 'String', cVal);
            end
            
%             set(eval(strV), 'Visible', 'on', 'String', handles.userdata.ScaleObjs.default{i});

            if (handles.userdata.ScaleObjs.args(i) >= 100)
                set(eval(strVC), 'Visible', 'on');
            else
                set(eval(strVC), 'Visible', 'off');
            end
        end
        for j = length(handles.userdata.ScaleObjs.vars)+1:10
            strVT = sprintf('handles.vt%d',j);
            strV = sprintf('handles.v%d', j);
            set(eval(strVT), 'Visible', 'off', 'String', '');
            set(eval(strV), 'Visible', 'off', 'String', '');
            strVC = sprintf('handles.vc%d', j);
            set(eval(strVC), 'Visible', 'off');
        end
        set(handles.IOT_tag, 'Position', handles.userdata.IOT_tag_pos+[0 (10-length(handles.userdata.ScaleObjs.vars))*height 0 0]);
        set(handles.IO_tag, 'Position', handles.userdata.IO_tag_pos+[0 (10-length(handles.userdata.ScaleObjs.vars))*height 0 0]);
        set(handles.IOL_tag, 'Position', handles.userdata.IOL_tag_pos+[0 0 0 (10-length(handles.userdata.ScaleObjs.vars))*height]);
        implemented = 1;
    case 'reg2procalign_meta'
        for i = 1:length(handles.userdata.reg2procalign_meta.vars)
            strVT = sprintf('handles.vt%d',i);
            strV = sprintf('handles.v%d', i);
            strVC = sprintf('handles.vc%d', i);
            set(eval(strVT), 'Visible', 'on', 'String', handles.userdata.reg2procalign_meta.vars{i});
            if strcmp(handles.userdata.reg2procalign_meta.vars{i},'SampleMesh');
                set(eval(strV), 'Visible', 'on', 'String', handles.userdata.mesh(2:end),'Style','popupmenu');
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.reg2procalign_meta.%s=contentStr{valsMenu};',handles.userdata.reg2procalign_meta.vars{i}));
            else
                cValStr = ['handles.userdata.reg2procalign_meta.' handles.userdata.reg2procalign_meta.vars{i}];
                cVal = eval(cValStr);
                if handles.userdata.reg2procalign_meta.args(i) < 10
                    set(eval(strV), 'Visible', 'on', 'String', num2str(cVal));
                elseif handles.userdata.reg2procalign_meta.args(i) >= 10
                    set(eval(strV), 'Visible', 'on', 'String', cVal);
                end
%                 set(eval(strV), 'Visible', 'on', 'String', handles.userdata.reg2procalign_meta.default{i});
            end

            if (handles.userdata.reg2procalign_meta.args(i) >= 100)
                set(eval(strVC), 'Visible', 'on');
            else
                set(eval(strVC), 'Visible', 'off');
            end
        end
        for j = length(handles.userdata.reg2procalign_meta.vars)+1:10
            strVT = sprintf('handles.vt%d',j);
            strV = sprintf('handles.v%d', j);
            set(eval(strVT), 'Visible', 'off', 'String', '');
            set(eval(strV), 'Visible', 'off', 'String', '');
            strVC = sprintf('handles.vc%d', j);
            set(eval(strVC), 'Visible', 'off');
        end
        set(handles.IOT_tag, 'Position', handles.userdata.IOT_tag_pos+[0 (10-length(handles.userdata.reg2procalign_meta.vars))*height 0 0]);
        set(handles.IO_tag, 'Position', handles.userdata.IO_tag_pos+[0 (10-length(handles.userdata.reg2procalign_meta.vars))*height 0 0]);
        set(handles.IOL_tag, 'Position', handles.userdata.IOL_tag_pos+[0 0 0 (10-length(handles.userdata.reg2procalign_meta.vars))*height]);
        implemented = 1;
    case 'InHouse_Fix'
        for i = 1:length(handles.userdata.InHouse_Fix.vars)
            strVT = sprintf('handles.vt%d',i);
            strV = sprintf('handles.v%d', i);
            strVC = sprintf('handles.vc%d', i);
            set(eval(strVT), 'Visible', 'on', 'String', handles.userdata.InHouse_Fix.vars{i});

            if strcmp(handles.userdata.InHouse_Fix.vars{i},'Mesh')
                set(eval(strV), 'Visible', 'on', 'String', handles.userdata.mesh,'Style','popupmenu');
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.InHouse_Fix.%s=contentStr{valsMenu};',handles.userdata.InHouse_Fix.vars{i}));
            elseif strcmp(handles.userdata.InHouse_Fix.vars{i},'Connectivity')
                set(eval(strV), 'Visible', 'on', 'String', {'(6+,18)';'(18,6+)';'(6,26)';'(26,6)'},'Style','popupmenu');
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.InHouse_Fix.%s=contentStr{valsMenu};',handles.userdata.InHouse_Fix.vars{i}));
            else
                cValStr = ['handles.userdata.InHouse_Fix.' handles.userdata.InHouse_Fix.vars{i}];
                cVal = eval(cValStr);
                if handles.userdata.InHouse_Fix.args(i) < 10
                    set(eval(strV), 'Visible', 'on', 'String', num2str(cVal));
                elseif handles.userdata.InHouse_Fix.args(i) >= 10
                    set(eval(strV), 'Visible', 'on', 'String', cVal);
                end                
%                set(eval(strV), 'Visible', 'on', 'String', handles.userdata.InHouse_Fix.default{i});
            end

            if (handles.userdata.InHouse_Fix.args(i) >= 100)
                set(eval(strVC), 'Visible', 'on');
            else
                set(eval(strVC), 'Visible', 'off');
            end
        end
        for j = length(handles.userdata.InHouse_Fix.vars)+1:10
            strVT = sprintf('handles.vt%d',j);
            strV = sprintf('handles.v%d', j);
            set(eval(strVT), 'Visible', 'off', 'String', '');
            set(eval(strV), 'Visible', 'off', 'String', '');
            strVC = sprintf('handles.vc%d', j);
            set(eval(strVC), 'Visible', 'off');
        end
        set(handles.IOT_tag, 'Position', handles.userdata.IOT_tag_pos+[0 (10-length(handles.userdata.InHouse_Fix.vars))*height 0 0]);
        set(handles.IO_tag, 'Position', handles.userdata.IO_tag_pos+[0 (10-length(handles.userdata.InHouse_Fix.vars))*height 0 0]);
        set(handles.IOL_tag, 'Position', handles.userdata.IOL_tag_pos+[0 0 0 (10-length(handles.userdata.InHouse_Fix.vars))*height]);
        implemented = 1;
    case 'PDM_Fix'
        for i = 1:length(handles.userdata.PDM_Fix.vars)
            strVT = sprintf('handles.vt%d',i);
            strV = sprintf('handles.v%d', i);
            strVC = sprintf('handles.vc%d', i);
            set(eval(strVT), 'Visible', 'on', 'String', handles.userdata.PDM_Fix.vars{i});

            if strcmp(handles.userdata.PDM_Fix.vars{i},'Mesh');
                set(eval(strV), 'Visible', 'on', 'String', handles.userdata.mesh,'Style','popupmenu');
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.PDM_Fix.%s=contentStr{valsMenu};',handles.userdata.PDM_Fix.vars{i}));
            else
                cValStr = ['handles.userdata.PDM_Fix.' handles.userdata.PDM_Fix.vars{i}];
                cVal = eval(cValStr);
                if handles.userdata.PDM_Fix.args(i) < 10
                    set(eval(strV), 'Visible', 'on', 'String', num2str(cVal));
                elseif handles.userdata.PDM_Fix.args(i) >= 10
                    set(eval(strV), 'Visible', 'on', 'String', cVal);
                end
                
%                 set(eval(strV), 'Visible', 'on', 'String', handles.userdata.PDM_Fix.default{i});
            end

            if (handles.userdata.PDM_Fix.args(i) >= 100)
                set(eval(strVC), 'Visible', 'on');
            else
                set(eval(strVC), 'Visible', 'off');
            end
        end
        for j = length(handles.userdata.PDM_Fix.vars)+1:10
            strVT = sprintf('handles.vt%d',j);
            strV = sprintf('handles.v%d', j);
            set(eval(strVT), 'Visible', 'off', 'String', '');
            set(eval(strV), 'Visible', 'off', 'String', '');
            strVC = sprintf('handles.vc%d', j);
            set(eval(strVC), 'Visible', 'off');
        end
        set(handles.IOT_tag, 'Position', handles.userdata.IOT_tag_pos+[0 (10-length(handles.userdata.PDM_Fix.vars))*height 0 0]);
        set(handles.IO_tag, 'Position', handles.userdata.IO_tag_pos+[0 (10-length(handles.userdata.PDM_Fix.vars))*height 0 0]);
        set(handles.IOL_tag, 'Position', handles.userdata.IOL_tag_pos+[0 0 0 (10-length(handles.userdata.PDM_Fix.vars))*height]);
        implemented = 1;
    case 'Config'
        for i = 1:length(handles.userdata.Config.vars)
            strVT = sprintf('handles.vt%d',i);
            strV = sprintf('handles.v%d', i);
            strVC = sprintf('handles.vc%d', i);
            set(eval(strVT), 'Visible', 'on', 'String', handles.userdata.Config.vars{i});
            defaultVal = sprintf('handles.userdata.Config.%s', handles.userdata.Config.vars{i});
            set(eval(strV), 'Visible', 'on', 'String', eval(defaultVal));
            if (handles.userdata.Config.args(i) >= 100)
                set(eval(strVC), 'Visible', 'on');
            else
                set(eval(strVC), 'Visible', 'off');
            end
        end
%         set(handles.vc1, 'Visible', 'on');
%         set(handles.vc2, 'Visible', 'on');
        for j = length(handles.userdata.Config.vars)+1:10
            strVT = sprintf('handles.vt%d',j);
            strV = sprintf('handles.v%d', j);
            set(eval(strVT), 'Visible', 'off', 'String', '');
            set(eval(strV), 'Visible', 'off', 'String', '');
            strVC = sprintf('handles.vc%d', j);
            set(eval(strVC), 'Visible', 'off');
        end
        set(handles.IOT_tag, 'Position', handles.userdata.IOT_tag_pos);
        set(handles.IO_tag, 'Position', handles.userdata.IO_tag_pos);
        set(handles.IOL_tag, 'Position', handles.userdata.IOL_tag_pos);
        implemented = 1;
    case 'Import'
        for i = 1:length(handles.userdata.Import.vars)
            strVT = sprintf('handles.vt%d',i);
            strV = sprintf('handles.v%d', i);
            strVC = sprintf('handles.vc%d', i);
            set(eval(strVT), 'Visible', 'on', 'String', handles.userdata.Import.vars{i});

            if strcmp(handles.userdata.Import.vars{i},'Mesh');
                set(eval(strV), 'Visible', 'on', 'String', handles.userdata.mesh,'Style','popupmenu');
                contentStr = get(eval(strV),'String');
                valsMenu = get(eval(strV),'Value');
                eval(sprintf('handles.userdata.Import.%s=contentStr{valsMenu};',handles.userdata.Import.vars{i}));
            else
                cValStr = ['handles.userdata.Import.' handles.userdata.Import.vars{i}];
                cVal = eval(cValStr);
                if handles.userdata.Import.args(i) < 10
                    set(eval(strV), 'Visible', 'on', 'String', num2str(cVal));
                elseif handles.userdata.Import.args(i) >= 10
                    set(eval(strV), 'Visible', 'on', 'String', cVal);
                end
%                 set(eval(strV), 'Visible', 'on', 'String', handles.userdata.Import.default{i});
            end

            if (handles.userdata.Import.args(i) >= 100)
                set(eval(strVC), 'Visible', 'on');
            else
                set(eval(strVC), 'Visible', 'off');
            end
        end
        for j = length(handles.userdata.Import.vars)+1:10
            strVT = sprintf('handles.vt%d',j);
            strV = sprintf('handles.v%d', j);
            set(eval(strVT), 'Visible', 'off', 'String', '');
            set(eval(strV), 'Visible', 'off', 'String', '');
            strVC = sprintf('handles.vc%d', j);
            set(eval(strVC), 'Visible', 'off');
        end
        set(handles.IOT_tag, 'Position', handles.userdata.IOT_tag_pos+[0 (10-length(handles.userdata.Import.vars))*height 0 0]);
        set(handles.IO_tag, 'Position', handles.userdata.IO_tag_pos+[0 (10-length(handles.userdata.Import.vars))*height 0 0]);
        set(handles.IOL_tag, 'Position', handles.userdata.IOL_tag_pos+[0 0 0 (10-length(handles.userdata.Import.vars))*height]);
        implemented = 1;
    case {'bim2gipl', 'fix2gipl', 'gipl2bim', 'smo2surf_para_meta', 'surf_para_meta2smo', 'des2meta_coef', 'meta_coef2des', 'ellalign_meta_coef2reg'}
        set(handles.vt1, 'Visible', 'on', 'String', 'OutDirectory');
        
        st1 = ['handles.userdata.' handles.userdata.keyword.selected]
        st2 = '.vars{1}';
        cValStr = [st1 st2];
        cVal = eval([st1 '.' eval(cValStr)]);
        if eval([st1 '.args(1)']) < 10
            set(handles.v1, 'Visible', 'on', 'String', num2str(cVal));
        elseif eval([st1 '.args(1)']) >= 10
            set(handles.v1, 'Visible', 'on', 'String', cVal);
        end
%        set(handles.v1, 'Visible', 'on', 'String', './');
        
        set(handles.vc1, 'Visible', 'on');
            
        for j = 2:10
            strVT = sprintf('handles.vt%d',j);
            strV = sprintf('handles.v%d', j);
            set(eval(strVT), 'Visible', 'off', 'String', '');
            set(eval(strV), 'Visible', 'off', 'String', '');
            strVC = sprintf('handles.vc%d', j);
            set(eval(strVC), 'Visible', 'off');
        end
        set(handles.IOT_tag, 'Position', handles.userdata.IOT_tag_pos+[0 9*height 0 0]);
        set(handles.IO_tag, 'Position', handles.userdata.IO_tag_pos+[0 9*height 0 0]);
        set(handles.IOL_tag, 'Position', handles.userdata.IOL_tag_pos+[0 0 0 9*height]);
        implemented = 1;
    otherwise
        for i = 1:10
            strVT = sprintf('handles.vt%d',i);
            strV = sprintf('handles.v%d', i);
            strVC = sprintf('handles.vc%d', i);
            set(eval(strVT), 'Visible', 'off');
            set(eval(strV), 'Visible', 'off', 'String', '');
            set(eval(strVC), 'Visible', 'off');
        end
        set(handles.IOT_tag, 'Position', handles.userdata.IOT_tag_pos+[0 10*height 0 0]);
        set(handles.IO_tag, 'Position', handles.userdata.IO_tag_pos+[0 10*height 0 0]);
        set(handles.IOL_tag, 'Position', handles.userdata.IOL_tag_pos+[0 0 0 10*height]);
        implemented = 1;
end

if ~strcmp(handles.userdata.keyword.selected, 'Config')
    set(handles.IOT_tag, 'Visible', 'on');
    set(handles.IO_tag, 'Visible', 'on');
    set(handles.IOL_tag, 'Visible', 'on', 'String', '');
else
    set(handles.IOT_tag, 'Visible', 'off');
    set(handles.IO_tag, 'Visible', 'off');
    set(handles.IOL_tag, 'Visible', 'off');
end

if (implemented == 0)
    set(handles.RunTag, 'Enable', 'off');
else
    set(handles.RunTag, 'Enable', 'on');
end

return;