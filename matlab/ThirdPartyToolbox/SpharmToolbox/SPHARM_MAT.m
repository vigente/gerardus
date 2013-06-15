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

function varargout = SPHARM_MAT(varargin)
% SPHARM_MAT M-file for SPHARM_MAT.fig
%      SPHARM_MAT, by itself, creates a new SPHARM_MAT or raises the existing
%      singleton*.
%
%      H = SPHARM_MAT returns the handle to a new SPHARM_MAT or the handle to
%      the existing singleton*.
%
%      SPHARM_MAT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPHARM_MAT.M with the given input arguments.
%
%      SPHARM_MAT('Property','Value',...) creates a new SPHARM_MAT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SPHARM_MAT_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SPHARM_MAT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help SPHARM_MAT

% Last Modified by GUIDE v2.5 16-Apr-2009 11:12:27

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SPHARM_MAT_OpeningFcn, ...
                   'gui_OutputFcn',  @SPHARM_MAT_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

return;

% --- Executes just before SPHARM_MAT is made visible.
function SPHARM_MAT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SPHARM_MAT (see VARARGIN)

% Choose default command line output for SPHARM_MAT
%global userdata;

% The below routing doesnot work
%initializeParams(hObject, handles);

% Merging the entire initialization steps into this!!
% % Start initialization

% List of input object names
handles.userdata.inObjs = {};
handles.userdata.inObjs2 = {};

handles.userdata = initOptions(handles.userdata);

% Default position of three GUI components
handles.userdata.IOT_tag_pos = get(handles.IOT_tag, 'Position');
handles.userdata.IO_tag_pos = get(handles.IO_tag, 'Position');
handles.userdata.IOL_tag_pos = get(handles.IOL_tag, 'Position');

set(handles.IOT2_tag, 'Position', handles.userdata.IOT_tag_pos);
set(handles.IO2_tag, 'Position', handles.userdata.IO_tag_pos);
set(handles.IOL2_tag, 'Position', handles.userdata.IOL_tag_pos);

handles.userdata.cpath = pwd;
if exist(fullfile(handles.userdata.cpath,'SPHARM_MAT.m'),'file')
    addpath(handles.userdata.cpath);
end

if (handles.userdata.Config.loaded == 0) | (handles.userdata.ExpAlig.exist == 0)
    set(handles.ExpNAliTag, 'Enable', 'off');
end
% End initialization

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SPHARM_MAT wait for user response (see UIRESUME)
% uiwait(handles.figure1);

%end

return;

% --- Outputs from this function are returned to the command line.
function varargout = SPHARM_MAT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%rmpath(handles.userdata.cpath);
varargout{1} = handles.output;

%end
return;


% --- Executes on button press in Quit.
function Quit_Callback(hObject, eventdata, handles)
% hObject    handle to Quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Added by Sungeun Kim _Sep_08_08
% Before close the current GUI window, all variables which are created in
% the current study should be clear.
% if handles.userdata.config.loaded == 1
%     rmpath(handles.userdata.config.PDM_path);
%     rmpath(handles.userdata.config.Lapack_path);
% end
delete(get(0,'Children'));
clc;
%end
return;


% --- Executes on button press in InitParam.
function InitParam_Callback(hObject, eventdata, handles)
% hObject    handle to InitParam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global userdata;

state = get(hObject, 'Value');
set(handles.tools_tag, 'Value', 1);
set(handles.utils_tag, 'Value', 1);

if (state == 1)
%    handles.userdata.selected = 'Parametrization';
    
    hs = findobj('Style', 'togglebutton', '-and', 'Value', 1);
    set(hs, 'ForegroundColor', [0, 0, 0], 'Value', 0);
    
    set(handles.InitParam, 'ForegroundColor', [1.0, 0, 0], 'Value', 1);
    set(handles.OptionTitle, 'String', 'Parameterization');
    set(handles.LoadTag, 'Visible', 'on');
    set(handles.SaveTag, 'Visible', 'on');
    set(handles.RunTag, 'Visible', 'on', 'String', 'Run');
    set(handles.CancelTag, 'Visible', 'on');

    set(handles.MethodText, 'Visible', 'on');
    set(handles.MethodPopup, 'Visible', 'on', 'String',{'CALD';'PDM'}, 'Value', 1);
    
    handles.userdata.keyword.selected = 'ParamCALD';
    eraseOptions(handles);
    handles = dispOptions(handles);
    disp(handles.userdata.keyword.selected);
else
%    handles.userdata.selected = '';
    handles.userdata.keyword.selected = '';
    eraseOptions(handles);
    
    set(handles.InitParam, 'ForegroundColor', [0.0, 0.0, 0.0]);    
    set(handles.OptionTitle, 'String', '');    
    set(handles.LoadTag, 'Visible', 'off');
    set(handles.SaveTag, 'Visible', 'off');
    set(handles.RunTag, 'Visible', 'off');
    set(handles.CancelTag, 'Visible', 'off');

    set(handles.MethodText, 'Visible', 'off');
    set(handles.MethodPopup, 'Visible', 'off');
end    

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes on button press in Expansion.
function Expansion_Callback(hObject, eventdata, handles)
% hObject    handle to Expansion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global userdata;

state = get(hObject, 'Value');
set(handles.tools_tag, 'Value', 1);
set(handles.utils_tag, 'Value', 1);

if (state == 1)
%    handles.userdata.selected = 'Expansion';
    
    hs = findobj('Style', 'togglebutton', '-and', 'Value', 1);
    set(hs, 'ForegroundColor', [0, 0, 0], 'Value', 0);
    
    set(handles.Expansion, 'ForegroundColor', [1.0, 0, 0], 'Value', 1);
    set(handles.OptionTitle, 'String', 'SPHARM Expansion');        
    set(handles.LoadTag, 'Visible', 'on');
    set(handles.SaveTag, 'Visible', 'on');
    set(handles.RunTag, 'Visible', 'on', 'String', 'Run');
    set(handles.CancelTag, 'Visible', 'on');

    set(handles.MethodText, 'Visible', 'on');
    set(handles.MethodPopup, 'Visible', 'on', 'String',{'LSF'}, 'Value', 1);
    
    handles.userdata.keyword.selected = 'ExpLSF';
    eraseOptions(handles);
    handles = dispOptions(handles);
    disp(handles.userdata.keyword.selected);
else
%    handles.userdata.selected = '';
    handles.userdata.keyword.selected = '';
    eraseOptions(handles);
    
    set(handles.Expansion, 'ForegroundColor', [0.0, 0.0, 0.0]);    
    set(handles.OptionTitle, 'String', '');        
    set(handles.LoadTag, 'Visible', 'off');
    set(handles.SaveTag, 'Visible', 'off');
    set(handles.RunTag, 'Visible', 'off');
    set(handles.CancelTag, 'Visible', 'off');

    set(handles.MethodText, 'Visible', 'off');
    set(handles.MethodPopup, 'Visible', 'off');

end    

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes on button press in Alignment.
function Alignment_Callback(hObject, eventdata, handles)
% hObject    handle to Alignment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Alignment
%global userdata;

state = get(hObject, 'Value');
set(handles.tools_tag, 'Value', 1);
set(handles.utils_tag, 'Value', 1);

if (state == 1)
%    handles.userdata.selected = 'Alignment';    
    
    hs = findobj('Style', 'togglebutton', '-and', 'Value', 1);
    set(hs, 'ForegroundColor', [0, 0, 0], 'Value', 0);
    
    set(handles.Alignment, 'ForegroundColor', [1.0, 0, 0], 'Value', 1);
    set(handles.OptionTitle, 'String', 'Alignment');        
    set(handles.LoadTag, 'Visible', 'on');
    set(handles.SaveTag, 'Visible', 'on');
    set(handles.RunTag, 'Visible', 'on', 'String', 'Run');
    set(handles.CancelTag, 'Visible', 'on');
    
    set(handles.MethodText, 'Visible', 'on');
    set(handles.MethodPopup, 'Visible', 'on', 'String',{'SHREC';'FOE'}, 'Value', 1);
    
    handles.userdata.keyword.selected = 'AligSHREC';
    eraseOptions(handles);
    handles = dispOptions(handles);
    disp(handles.userdata.keyword.selected);
else
%    handles.userdata.selected = '';
    handles.userdata.keyword.selected = '';
    eraseOptions(handles);
    
    set(handles.Alignment, 'ForegroundColor', [0.0, 0.0, 0.0]);    
    set(handles.OptionTitle, 'String', '');        
    set(handles.LoadTag, 'Visible', 'off');
    set(handles.SaveTag, 'Visible', 'off');
    set(handles.RunTag, 'Visible', 'off');
    set(handles.CancelTag, 'Visible', 'off');

    set(handles.MethodText, 'Visible', 'off');
    set(handles.MethodPopup, 'Visible', 'off');

end    

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes on button press in ExpNAliTag.
function ExpNAliTag_Callback(hObject, eventdata, handles)
% hObject    handle to ExpNAliTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ExpNAliTag
%global userdata;

state = get(hObject, 'Value');
set(handles.tools_tag, 'Value', 1);
set(handles.utils_tag, 'Value', 1);

if (state == 1)
%    handles.userdata.keyword.selected = 'ExpNAli';  

    hs = findobj('Style', 'togglebutton', '-and', 'Value', 1);
    set(hs, 'ForegroundColor', [0, 0, 0], 'Value', 0);

    set(handles.MethodText, 'Visible', 'off');
    set(handles.MethodPopup, 'Visible', 'off');    
    
    set(handles.ExpNAliTag, 'ForegroundColor', [1.0, 0, 0], 'Value', 1);
    set(handles.OptionTitle, 'String', 'Expansion & Alignment (SHPARM-PDM)');          

    set(handles.LoadTag, 'Visible', 'on');
    set(handles.SaveTag, 'Visible', 'on');
    set(handles.RunTag, 'Visible', 'on', 'String', 'OK');
    set(handles.CancelTag, 'Visible', 'on');

    handles.userdata.keyword.selected = 'ExpAlig';
    eraseOptions(handles);
    handles = dispOptions(handles);
    disp(handles.userdata.keyword.selected);    
else
%    handles.userdata.selected = '';    
    handles.userdata.keyword.selected = '';
    eraseOptions(handles);
    
    set(handles.ExpNAliTag, 'ForegroundColor', [0.0, 0.0, 0.0]);    
    set(handles.OptionTitle, 'String', '');          
    set(handles.LoadTag, 'Visible', 'off');
    set(handles.SaveTag, 'Visible', 'off');
    set(handles.RunTag, 'Visible', 'off', 'String', 'Run');
    set(handles.CancelTag, 'Visible', 'off');
end    

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes on button press in StatAnalysis.
function StatAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to StatAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of StatAnalysis
%global userdata;

state = get(hObject, 'Value');
set(handles.tools_tag, 'Value', 1);
set(handles.utils_tag, 'Value', 1);

if (state == 1)
%    handles.userdata.selected = 'StatAnalysis';    
    
    hs = findobj('Style', 'togglebutton', '-and', 'Value', 1);
    set(hs, 'ForegroundColor', [0, 0, 0], 'Value', 0);
    
    set(handles.StatAnalysis, 'ForegroundColor', [1.0, 0, 0], 'Value', 1);
    set(handles.OptionTitle, 'String', 'Statistical Analysis');          
    set(handles.LoadTag, 'Visible', 'on');
    set(handles.SaveTag, 'Visible', 'on');
    set(handles.RunTag, 'Visible', 'on', 'String', 'Run');
    set(handles.CancelTag, 'Visible', 'on');

    set(handles.MethodText, 'Visible', 'on');
    set(handles.MethodPopup, 'Visible', 'on', 'String',{'t_map';'PCA'}, 'Value', 1);

    handles.userdata.keyword.selected = 't_map';
    eraseOptions(handles);
    handles = dispOptions(handles);
    disp(handles.userdata.keyword.selected);
else
%    handles.userdata.selected = '';
    handles.userdata.keyword.selected = '';
    eraseOptions(handles);
    
    set(handles.StatAnalysis, 'ForegroundColor', [0.0, 0.0, 0.0]);    
    set(handles.OptionTitle, 'String', '');          
    set(handles.LoadTag, 'Visible', 'off');
    set(handles.SaveTag, 'Visible', 'off');
    set(handles.RunTag, 'Visible', 'off');
    set(handles.CancelTag, 'Visible', 'off');
    
    set(handles.MethodText, 'Visible', 'off');
    set(handles.MethodPopup, 'Visible', 'off');
    
end    

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes on button press in Import_tag.
function Import_tag_Callback(hObject, eventdata, handles)
% hObject    handle to Import_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Import_tag

state = get(hObject, 'Value');
set(handles.tools_tag, 'Value', 1);
set(handles.utils_tag, 'Value', 1);

if (state == 1)
%    handles.userdata.keyword.selected = 'Import';    

    hs = findobj('Style', 'togglebutton', '-and', 'Value', 1);
    set(hs, 'ForegroundColor', [0, 0, 0], 'Value', 0);

    set(handles.LoadTag, 'Visible', 'off');
    set(handles.SaveTag, 'Visible', 'off');    
    
    set(handles.MethodText, 'Visible', 'off');
    set(handles.MethodPopup, 'Visible', 'off');
    
    set(handles.Import_tag, 'ForegroundColor', [1.0, 0, 0], 'Value', 1);
    set(handles.OptionTitle, 'String', 'Import object volumes');          
    set(handles.RunTag, 'Visible', 'on', 'String', 'OK');
    set(handles.CancelTag, 'Visible', 'on');
    
    handles.userdata.keyword.selected = 'Import';
    eraseOptions(handles);
    handles = dispOptions(handles);
    disp(handles.userdata.keyword.selected);    
else
%    handles.userdata.selected = '';    
    handles.userdata.keyword.selected = '';
    eraseOptions(handles);
    
    set(handles.Import_tag, 'ForegroundColor', [0.0, 0.0, 0.0]);    
    set(handles.OptionTitle, 'String', '');          
    set(handles.RunTag, 'Visible', 'off', 'String', 'Run');
    set(handles.CancelTag, 'Visible', 'off');
end    

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes on button press in DisplayTag.
function DisplayTag_Callback(hObject, eventdata, handles)
% hObject    handle to DisplayTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of DisplayTag
%global userdata;

state = get(hObject, 'Value');
set(handles.tools_tag, 'Value', 1);
set(handles.utils_tag, 'Value', 1);

if (state == 1)
%    handles.userdata.keyword.selected = 'DisplayRes';    
    
    hs = findobj('Style', 'togglebutton', '-and', 'Value', 1);
    set(hs, 'ForegroundColor', [0, 0, 0], 'Value', 0);
    
    set(handles.LoadTag, 'Visible', 'off');
    set(handles.SaveTag, 'Visible', 'off');    

    set(handles.MethodText, 'Visible', 'on');
    set(handles.MethodPopup, 'Visible', 'on', 'String',{'res_t_map';'res_PCA'}, 'Value', 1);
    
    set(handles.DisplayTag, 'ForegroundColor', [1.0, 0, 0], 'Value', 1);
    set(handles.OptionTitle, 'String', 'Display Stat Results');          
    set(handles.RunTag, 'Visible', 'on', 'String', 'OK');
    set(handles.CancelTag, 'Visible', 'on');
    
    handles.userdata.keyword.selected = 'res_t_map';
    eraseOptions(handles);
    handles = dispOptions(handles);
    disp(handles.userdata.keyword.selected);    
else
%    handles.userdata.selected = '';
    handles.userdata.keyword.selected = '';
    eraseOptions(handles);
    
    set(handles.DisplayTag, 'ForegroundColor', [0.0, 0.0, 0.0]);    
    set(handles.OptionTitle, 'String', '');          
    set(handles.RunTag, 'Visible', 'off', 'String', 'Run');
    set(handles.CancelTag, 'Visible', 'off');
    
    set(handles.MethodText, 'Visible', 'off');
    set(handles.MethodPopup, 'Visible', 'off');
    
end    

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes on button press in LoadTag.
function LoadTag_Callback(hObject, eventdata, handles)
% hObject    handle to LoadTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global userdata;

[inFile, inPath, fIndex] = uigetfile({'*.txt','Text Config Files (*.txt)'; ... 
            '*.*',  'All Files (*.*)'}, ...
            'Select a configuration file');

if fIndex ~=0        
    cmdStr1 = sprintf('handles.userdata.%s=loadOptions([inPath, inFile], handles.userdata.%s);',handles.userdata.keyword.selected, handles.userdata.keyword.selected);
    cmdStr2 = sprintf('updateOptions(handles, handles.userdata.%s);',handles.userdata.keyword.selected);

    eval(cmdStr1);
    [handles.userdata.inObjs, handles.userdata.inObjs2] = LoadObjList([inPath,inFile]);
    eval(cmdStr2);
end

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes on button press in SaveTag.
function SaveTag_Callback(hObject, eventdata, handles)
% hObject    handle to SaveTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global userdata;

[outFile, outPath, fIndex] = uiputfile({'*.txt','Text Config Files (*.txt)'; ... 
            '*.*',  'All Files (*.*)'}, ...
            'Save Configuration As');

if fIndex ~= 0        

    cmdStr1 = sprintf('saveOptions([outPath, outFile], handles.userdata.%s, handles.userdata.inObjs, handles.userdata.inObjs2);',handles.userdata.keyword.selected);
    eval(cmdStr1);

end

return;


% --- Executes on button press in RunTag.
function RunTag_Callback(hObject, eventdata, handles)
% hObject    handle to RunTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global userdata;

disp(sprintf('%s', handles.userdata.keyword.selected));

%updateVariables(hObject, eventdata, handles);

switch handles.userdata.keyword.selected
    case 'ParamCALD'
        SpharmMatParameterization(handles.userdata.ParamCALD, handles.userdata.inObjs, 'ParamCALD');
    case 'ParamPDM'
        if (handles.userdata.Config.loaded == 0)
            errordlg('Path to PDM program is not set.\n You need to set the path to PDM program in config menu of tool popup menu.');
        elseif (handles.userdata.ParamPDM.exist == 0)
            errordlg(sprintf('%s function does not exist in the PDM directory, %s',handles.userdata.ParamPDM.command, handles.userdata.Config.PDM_path));
        else
            SpharmMatParameterization(handles.userdata.ParamPDM, handles.userdata.inObjs, 'ParamPDM');
        end
    case 'ParamQuad'
        SpharmMatParameterization(handles.userdata.ParamQuad, handles.userdata.inObjs, 'ParamQuad');        
    case 'ParamConf'
        SpharmMatParameterization(handles.userdata.ParamConf, handles.userdata.inObjs, 'ParamConf');        
    case 'ExpLSF'
        SpharmMatExpansion(handles.userdata.ExpLSF, handles.userdata.inObjs, 'ExpLSF');
    case 'ExpIRF'
        SpharmMatExpansion(handles.userdata.ExpIRF, handles.userdata.inObjs, 'ExpIRF');        
    case 'AligSHREC'
        SpharmMatAlignment(handles.userdata.AligSHREC, handles.userdata.inObjs, 'AligSHREC');
    case 'AligFOE'
        SpharmMatAlignment(handles.userdata.AligFOE, handles.userdata.inObjs, 'AligFOE');        
    case 'AligLandmark'
        SpharmMatAlignment(handles.userdata.AligLandmark, handles.userdata.inObjs, 'AligLandmark');        
    case 'ExpAlig'
        if (handles.userdata.Config.loaded == 0)
            errordlg('Path to PDM program is not set.\n You need to set the path to PDM program in config menu of tool popup menu.');
        elseif (handles.userdata.ExpAlig.exist == 0)
            errordlg(sprintf('%s function does not exist in the PDM directory, %s',handles.userdata.ExpAlig.command, handles.userdata.Config.PDM_path));
        else
            SpharmMatExpandAlignByPDM(handles.userdata.ExpAlig, handles.userdata.inObjs);
        end
    case 't_map'
        SpharmMatStatAnalysis(handles.userdata.t_map, handles.userdata.inObjs, handles.userdata.inObjs2, 't_map', handles.userdata.cpath);
    case 'PCA'
        SpharmMatStatAnalysis(handles.userdata.PCA, handles.userdata.inObjs, {}, 'PCA', handles.userdata.cpath);        
    case 'Import'
        if ~isempty(handles.userdata.inObjs)
            SpharmMatUtilImportObjs(handles.userdata.Import, handles.userdata.inObjs);
        end
    case 'bim2gipl'
        SpharmMatUtilFormatConvert(handles.userdata.bim2gipl, handles.userdata.inObjs, 'bim2gipl');
    case 'fix2gipl'
        SpharmMatUtilFormatConvert(handles.userdata.fix2gipl, handles.userdata.inObjs, 'fix2gipl');
    case 'gipl2bim'
        SpharmMatUtilFormatConvert(handles.userdata.gipl2bim, handles.userdata.inObjs, 'gipl2bim');
    case 'smo2surf_para_meta'
        SpharmMatUtilFormatConvert(handles.userdata.smo2surf_para_meta, handles.userdata.inObjs, 'smo2surf_para_meta');
    case 'surf_para_meta2smo'
        SpharmMatUtilFormatConvert(handles.userdata.surf_para_meta2smo, handles.userdata.inObjs, 'surf_para_meta2smo');
    case 'des2meta_coef'
        SpharmMatUtilFormatConvert(handles.userdata.des2meta_coef, handles.userdata.inObjs, 'des2meta_coef');
    case 'meta_coef2des'
        SpharmMatUtilFormatConvert(handles.userdata.meta_coef2des, handles.userdata.inObjs, 'meta_coef2des');
    case 'ellalign_meta_coef2reg'
        SpharmMatUtilFormatConvert(handles.userdata.ellalign_meta_coef2reg, handles.userdata.inObjs, 'ellalign_meta_coef2reg');        
    case 'reg2procalign_meta'
        SpharmMatUtilFormatConvert(handles.userdata.reg2procalign_meta, handles.userdata.inObjs, 'reg2procalign_meta');
    case 'InHouse_Fix'
        SpharmMatUtilTopologyFix(handles.userdata.InHouse_Fix, handles.userdata.inObjs, 'InHouse_Fix');
    case 'PDM_Fix'
        if (handles.userdata.Config.loaded == 0)
            errordlg('Path to PDM program is not set.\n You need to set the path to PDM program in config menu of tool popup menu.');
        elseif (handles.userdata.PDM_Fix.exist == 0)
            errordlg(sprintf('%s function does not exist in the PDM directory, %s',handles.userdata.PDM_Fix.command, handles.userdata.Config.PDM_path));
        else
            SpharmMatUtilTopologyFix(handles.userdata.PDM_Fix, handles.userdata.inObjs, 'PDM_Fix');
        end
    case 'res_t_map'
        h = figure('Name', 'Render Statistical Results','NumberTitle', 'off');
        cameratoolbar(h, 'Show');
        SpharmMatDisplayStat(handles.userdata.res_t_map, handles.userdata.inObjs, 'res_t_map', h, handles.userdata.cpath);
    case 'res_PCA'
        h = figure('Name', 'Render Statistical Results','NumberTitle', 'off');
        cameratoolbar(h, 'Show');
        SpharmMatDisplayStat(handles.userdata.res_PCA, handles.userdata.inObjs, 'res_PCA', h, handles.userdata.cpath);        
    case 'Config'
        handles.userdata.Config.PDM_path_pc = get(handles.v1, 'String');
        handles.userdata.Config.PDM_path_unix = get(handles.v2, 'String');        
        if ispc
            handles.userdata.Config.PDM_path = handles.userdata.Config.PDM_path_pc;
        else isunix
            handles.userdata.Config.PDM_path = handles.userdata.Config.PDM_path_unix;                        
        end
        handles.userdata.Config.loaded = 1;
        PDM_path_pc = handles.userdata.Config.PDM_path_pc;
        PDM_path_unix = handles.userdata.Config.PDM_path_unix;
        PDM_path = handles.userdata.Config.PDM_path;
        
        handles.userdata.ParamPDM.path = PDM_path;
        handles.userdata.PDM_Fix.path = PDM_path;
        handles.userdata.ExpAlig.path = PDM_path;
        handles.userdata.PDM_stat.path = PDM_path;
        guidata(hObject, handles);
        
        com1 = [PDM_path '/' handles.userdata.ParamPDM.command];
        com2 = [PDM_path '/' handles.userdata.ExpAlig.command];    
        com4 = [PDM_path '/' handles.userdata.PDM_Fix.command]; 

        if ispc
            com1 = [PDM_path '/' handles.userdata.ParamPDM.command '.exe'];
            com2 = [PDM_path '/' handles.userdata.ExpAlig.command '.exe'];    
            com4 = [PDM_path '/' handles.userdata.PDM_Fix.command '.exe']; 
        end
        if exist(com1, 'file')
            handles.userdata.ParamPDM.exist =1;
        else
            handles.userdata.ParamPDM.exist =0;
        end
        if exist(com2, 'file')
            handles.userdata.ExpAlig.exist =1;
        else
            handles.userdata.ExpAlig.exist =0;
        end
        if exist(com4, 'file')
            handles.userdata.PDM_Fix.exist =1;
        else
            handles.userdata.PDM_Fix.exist =0;
        end
        
        if (handles.userdata.Config.loaded == 1) & (handles.userdata.ExpAlig.exist == 1)
            set(handles.ExpNAliTag, 'Enable', 'on');
        else
            set(handles.ExpNAliTag, 'Enable', 'off');
        end
        
        save(handles.userdata.Config.filename, 'PDM_path_pc','PDM_path_unix');
    case 'DisplayObjs'
        
        SpharmMatUtilDisplayObjs(handles.userdata.DisplayObjs, handles.userdata.inObjs, handles.userdata.cpath);
    case 'AverageObjs'
        SpharmMatUtilAverageObjs(handles.userdata.AverageObjs, handles.userdata.inObjs)
    case 'ScaleObjs'
        SpharmMatUtilScaleObjs(handles.userdata.ScaleObjs, handles.userdata.inObjs);
end

CancelTag_Callback(hObject, eventdata, handles);

return;


% --- Executes on button press in CancelTag.
function CancelTag_Callback(hObject, eventdata, handles)
% hObject    handle to CancelTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global userdata;


hs = findobj('Style', 'togglebutton', '-and', 'Value', 1);
set(hs, 'ForegroundColor', [0, 0, 0], 'Value', 0);

set(get(handles.OptionTag, 'Children'), 'Visible', 'off');
set(handles.OptionTitle,'String', '', 'Visible', 'on');

set(handles.utils_tag, 'Value', 1);        
set(handles.tools_tag, 'Value', 1);
%handles.userdata.selected = '';
handles.userdata.keyword.selected = '';
eraseOptions(handles);

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes on selection change in MethodPopup.
function MethodPopup_Callback(hObject, eventdata, handles)
% hObject    handle to MethodPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns MethodPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from MethodPopup
%global userdata;
contents = get(hObject,'String');
selected = contents{get(hObject,'Value')};
eraseOptions(handles);

switch deblank(char(selected))
    case 'CALD'
        handles.userdata.keyword.selected = 'ParamCALD';
        disp(handles.userdata.keyword.selected);
    case 'PDM'
        handles.userdata.keyword.selected = 'ParamPDM';
        disp(handles.userdata.keyword.selected);   
    case 'InHouse_Quad'
        handles.userdata.keyword.selected = 'ParamQuad';
        disp(handles.userdata.keyword.selected);   
    case 'LSF'
        handles.userdata.keyword.selected = 'ExpLSF';
        disp(handles.userdata.keyword.selected);                
    case 'SHREC'
        handles.userdata.keyword.selected = 'AligSHREC';
        disp(handles.userdata.keyword.selected);                
    case 'FOE'
        handles.userdata.keyword.selected = 'AligFOE';
        disp(handles.userdata.keyword.selected);                
    otherwise
        handles.userdata.keyword.selected = deblank(char(selected));
        disp(handles.userdata.keyword.selected); 
end    

handles = dispOptions(handles);

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes during object creation, after setting all properties.
function MethodPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MethodPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


return;


% --- Executes on selection change in utils_tag.
function utils_tag_Callback(hObject, eventdata, handles)
% hObject    handle to utils_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns utils_tag contents as cell array
%        contents{get(hObject,'Value')} returns selected item from utils_tag

contents = get(hObject, 'String');
utils = contents(get(hObject, 'Value'));
set(handles.tools_tag, 'Value', 1);

hs = findobj('Style', 'togglebutton', '-and', 'Value', 1);
set(hs, 'ForegroundColor', [0, 0, 0], 'Value', 0);

set(handles.LoadTag, 'Visible', 'off');
set(handles.SaveTag, 'Visible', 'off');    

switch upper(deblank(char(utils)))
    case 'FORMATCONVERT'
        Convert_tag_Callback(hObject, eventdata, handles);
    case 'TOPOLOGYFIX'
        Topology_tag_Callback(hObject, eventdata, handles);
    case 'DISPLAYOBJS'
        Display_Obj_Callback(hObject, eventdata, handles);
    case 'AVERAGEOBJS'
        Average_Obj_Callback(hObject, eventdata, handles);
    case 'SCALEOBJS'
        Scale_Obj_Callback(hObject, eventdata, handles);
end
return;


% --- Executes during object creation, after setting all properties.
function utils_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to utils_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

return;


% --- Executes on button press in DisplayObjs.
function Display_Obj_Callback(hObject, eventdata, handles)
% hObject    handle to Convert_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Convert_tag

% hs = findobj('Style', 'togglebutton', '-and', 'Value', 1);
% set(hs, 'ForegroundColor', [0, 0, 0], 'Value', 0);

cSTR = {'DisplayObjs'};

set(handles.MethodText, 'Visible', 'off');
set(handles.MethodPopup, 'Visible', 'off');
    
set(handles.OptionTitle, 'String', 'Display Object');          
set(handles.RunTag, 'Visible', 'on', 'String', 'OK');
set(handles.CancelTag, 'Visible', 'on');

handles.userdata.keyword.selected = 'DisplayObjs';
eraseOptions(handles);
handles = dispOptions(handles);
disp(handles.userdata.keyword.selected);    

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes on button press in AverageObjs.
function Average_Obj_Callback(hObject, eventdata, handles)
% hObject    handle to Convert_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Convert_tag

% hs = findobj('Style', 'togglebutton', '-and', 'Value', 1);
% set(hs, 'ForegroundColor', [0, 0, 0], 'Value', 0);

cSTR = {'AverageObjs'};

set(handles.MethodText, 'Visible', 'off');
set(handles.MethodPopup, 'Visible', 'off', 'String',cSTR, 'Value', 1);

%    set(handles.Convert_tag, 'ForegroundColor', [1.0, 0, 0], 'Value', 1);
set(handles.OptionTitle, 'String', 'Create Average Object');          
set(handles.RunTag, 'Visible', 'on', 'String', 'OK');
set(handles.CancelTag, 'Visible', 'on');

handles.userdata.keyword.selected = 'AverageObjs';
eraseOptions(handles);
handles = dispOptions(handles);
disp(handles.userdata.keyword.selected);

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes on button press in ScaleObjs.
function Scale_Obj_Callback(hObject, eventdata, handles)
% hObject    handle to Convert_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Convert_tag

% hs = findobj('Style', 'togglebutton', '-and', 'Value', 1);
% set(hs, 'ForegroundColor', [0, 0, 0], 'Value', 0);

cSTR = {'ScaleObjs'};

set(handles.MethodText, 'Visible', 'off');
set(handles.MethodPopup, 'Visible', 'off', 'String',cSTR, 'Value', 1);

%    set(handles.Convert_tag, 'ForegroundColor', [1.0, 0, 0], 'Value', 1);
set(handles.OptionTitle, 'String', 'Scale Objects');          
set(handles.RunTag, 'Visible', 'on', 'String', 'OK');
set(handles.CancelTag, 'Visible', 'on');

handles.userdata.keyword.selected = 'ScaleObjs';
eraseOptions(handles);
handles = dispOptions(handles);
disp(handles.userdata.keyword.selected);

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes on button press in Convert_tag.
function Convert_tag_Callback(hObject, eventdata, handles)
% hObject    handle to Convert_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Convert_tag

% hs = findobj('Style', 'togglebutton', '-and', 'Value', 1);
% set(hs, 'ForegroundColor', [0, 0, 0], 'Value', 0);

cSTR = {'bim2gipl';'fix2gipl';'gipl2bim';'smo2surf_para_meta';'surf_para_meta2smo'; ...
    'des2meta_coef';'meta_coef2des';'ellalign_meta_coef2reg';'reg2procalign_meta'};

set(handles.MethodText, 'Visible', 'on');
set(handles.MethodPopup, 'Visible', 'on', 'String',cSTR, 'Value', 1);

%    set(handles.Convert_tag, 'ForegroundColor', [1.0, 0, 0], 'Value', 1);
set(handles.OptionTitle, 'String', 'File format conversion');          
set(handles.RunTag, 'Visible', 'on', 'String', 'OK');
set(handles.CancelTag, 'Visible', 'on');

handles.userdata.keyword.selected = 'bim2gipl';
eraseOptions(handles);
handles = dispOptions(handles);
disp(handles.userdata.keyword.selected);

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes on button press in Convert_tag.
function Topology_tag_Callback(hObject, eventdata, handles)
% hObject    handle to Convert_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Convert_tag

set(handles.MethodText, 'Visible', 'on');
set(handles.MethodPopup, 'Visible', 'on', 'String',{'InHouse_Fix';'PDM_Fix'}, 'Value', 1);

set(handles.LoadTag, 'Visible', 'on');
set(handles.SaveTag, 'Visible', 'on');    

%    set(handles.Convert_tag, 'ForegroundColor', [1.0, 0, 0], 'Value', 1);
set(handles.OptionTitle, 'String', 'Fix bad topology of binary objects');          
set(handles.RunTag, 'Visible', 'on', 'String', 'OK');
set(handles.CancelTag, 'Visible', 'on');


handles.userdata.keyword.selected = 'InHouse_Fix';
eraseOptions(handles);
handles = dispOptions(handles);
disp(handles.userdata.keyword.selected);

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes on selection change in tools_tag.
function tools_tag_Callback(hObject, eventdata, handles)
% hObject    handle to tools_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns tools_tag contents as cell array
%        contents{get(hObject,'Value')} returns selected item from tools_tag

contents = get(hObject, 'String');
tools = contents(get(hObject, 'Value'));
set(handles.utils_tag, 'Value', 1);

hs = findobj('Style', 'togglebutton', '-and', 'Value', 1);
set(hs, 'ForegroundColor', [0, 0, 0], 'Value', 0);

switch upper(deblank(char(tools)))
    case 'PWD'
        PWDTag_Callback(handles);
        set(hObject,'Value', 1);
    case 'CD'
        CDTag_Callback(handles);
        set(hObject,'Value', 1);        
    case 'CONFIG'
        Config_Callback(hObject, eventdata, handles);
end

return;


function PWDTag_Callback(handles)
% hObject    handle to PWDTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
eraseOptions(handles);

set(handles.LoadTag, 'Visible', 'off');
set(handles.SaveTag, 'Visible', 'off');    

set(handles.MethodText, 'Visible', 'off');
set(handles.MethodPopup, 'Visible', 'off');

set(handles.OptionTitle,'String', '');          
set(handles.RunTag, 'Visible', 'off');
set(handles.CancelTag, 'Visible', 'off');

cwd = sprintf('Current working directory: %s',pwd);
helpdlg(cwd);

return;


function CDTag_Callback(handles)
% hObject    handle to CDTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
eraseOptions(handles);

set(handles.LoadTag, 'Visible', 'off');
set(handles.SaveTag, 'Visible', 'off');    

set(handles.MethodText, 'Visible', 'off');
set(handles.MethodPopup, 'Visible', 'off');

set(handles.OptionTitle,'String', '');          
set(handles.RunTag, 'Visible', 'off');
set(handles.CancelTag, 'Visible', 'off');

new_dir = uigetdir(pwd);

if (new_dir ~= 0)
    cd(new_dir);

    cwd = sprintf('New working directory: %s',new_dir);
    helpdlg(cwd);
end

return;


% --- Executes on button press in ConfTag.
function Config_Callback(hObject, eventdata, handles)
% hObject    handle to ConfTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of ConfTag
% global userdata;

hs = findobj('Style', 'togglebutton', '-and', 'Value', 1);
set(hs, 'ForegroundColor', [0, 0, 0], 'Value', 0);

set(handles.LoadTag, 'Visible', 'off');
set(handles.SaveTag, 'Visible', 'off');    

set(handles.MethodText, 'Visible', 'off');
set(handles.MethodPopup, 'Visible', 'off');

%set(handles.ConfTag, 'ForegroundColor', [1.0, 0, 0], 'Value', 1);
set(handles.OptionTitle, 'String', 'Configuration');          
set(handles.RunTag, 'Visible', 'on', 'String', 'OK');
set(handles.CancelTag, 'Visible', 'on');

handles.userdata.keyword.selected = 'Config';
eraseOptions(handles);
handles = dispOptions(handles);
disp(handles.userdata.keyword.selected);

% Update handles structure
guidata(hObject, handles);

return;


function IOL_tag_Callback(hObject, eventdata, handles)
% hObject    handle to IOL_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IOL_tag as text
%        str2double(get(hObject,'String')) returns contents of IOL_tag as a double

handles.userdata.inObjs = get(hObject, 'String');

% Update handles structure
guidata(hObject, handles);
return;


% --- Executes during object creation, after setting all properties.
function IOL_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IOL_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

return;


% --- Executes on button press in IO_tag.
function IO_tag_Callback(hObject, eventdata, handles)
% hObject    handle to IO_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch handles.userdata.keyword.selected
    case 'ParamCALD'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.ParamCALD.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'ParamPDM'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.ParamPDM.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'ParamQuad'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.ParamQuad.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'ParamConf'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.ParamConf.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'ExpLSF'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.ExpLSF.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'ExpIRF'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.ExpIRF.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'AligSHREC'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.AligSHREC.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'AligFOE'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.AligFOE.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'AligLandmark'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.AligLandmark.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'ExpAlig'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.ExpAlig.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'PCA'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.PCA.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 't_map'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.t_map.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'Import'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.Import.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'res_t_map'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.res_t_map.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'off');
    case 'res_PCA'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.res_PCA.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'off');
    case 'bim2gipl'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.bim2gipl.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'fix2gipl'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.fix2gipl.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'gipl2bim'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.gipl2bim.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'smo2surf_para_meta'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.smo2surf_para_meta.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'surf_para_meta2smo'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.surf_para_meta2smo.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'des2meta_coef'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.des2meta_coef.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'meta_coef2des'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.meta_coef2des.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'ellalign_meta_coef2reg'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.ellalign_meta_coef2reg.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'reg2procalign_meta'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.reg2procalign_meta.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'InHouse_Fix'
%        setappdata(0, 'UseNativeSystemDialogs', 1);
        [filename, pathname, filterindex] = uigetfile(handles.userdata.InHouse_Fix.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'PDM_Fix'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.PDM_Fix.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'DisplayObjs'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.DisplayObjs.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'AverageObjs'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.AverageObjs.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'ScaleObjs'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.ScaleObjs.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    otherwise
        [filename, pathname, filterindex] = uigetfile(handles.userdata.others.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
end

if filterindex ~= 0
    if (class(filename) == 'cell') 
        d = size(filename,2);
        for i = 1:d
            files{i} = [pathname filename{i}];
        end
        set(handles.IOL_tag, 'String', files);

    else
        d = size(filename,1);
        for i = 1:d
            files{i} = [pathname filename(i,:)];
        end
        set(handles.IOL_tag, 'String', files);        
    end
    handles.userdata.inObjs = files;
end
% Update handles structure
guidata(hObject, handles);

return;


function v1_Callback(hObject, eventdata, handles)
% hObject    handle to v1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of v1 as text
%        str2double(get(hObject,'String')) returns contents of v1 as a double
st = get(handles.v1, 'Style');

if strcmp(st,'edit')
    vals = strtrim(get(handles.v1, 'String'));
elseif strcmp(st,'popupmenu')
    contents = get(handles.v1,'String');
    vals = strtrim(char(contents(get(handles.v1, 'Value'))));
end

Keyword1 = handles.userdata.keyword.selected;
varName = sprintf('handles.userdata.%s.vars{1}', Keyword1);
Keyword2 = eval(varName);
vars = sprintf('handles.userdata.%s.%s', Keyword1, Keyword2);

type = eval(sprintf('handles.userdata.%s.args(1)', Keyword1));

if type <10
    strC = sprintf('%s=str2num(vals);',vars);
    eval(strC);
else
    strC = sprintf('%s=vals;',vars);
    eval(strC);
end

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes during object creation, after setting all properties.
function v1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

return;


function v2_Callback(hObject, eventdata, handles)
% hObject    handle to v2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of v2 as text
%        str2double(get(hObject,'String')) returns contents of v2 as a double
st = get(handles.v2, 'Style');

if strcmp(st,'edit')
    vals = strtrim(get(handles.v2, 'String'));
elseif strcmp(st,'popupmenu')
    contents = get(handles.v2,'String');
    vals = strtrim(char(contents(get(handles.v2, 'Value'))));
end

Keyword1 = handles.userdata.keyword.selected;
varName = sprintf('handles.userdata.%s.vars{2}', Keyword1);
Keyword2 = eval(varName);
vars = sprintf('handles.userdata.%s.%s', Keyword1, Keyword2);

type = eval(sprintf('handles.userdata.%s.args(2)', Keyword1));

if type <10
    strC = sprintf('%s=str2num(vals);',vars);
    eval(strC);
else
    strC = sprintf('%s=vals;',vars);
    eval(strC);
end
    
% Update handles structure
guidata(hObject, handles);
return;


% --- Executes during object creation, after setting all properties.
function v2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


return;


function v3_Callback(hObject, eventdata, handles)
% hObject    handle to v3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of v3 as text
%        str2double(get(hObject,'String')) returns contents of v3 as a double
st = get(handles.v3, 'Style');

if strcmp(st,'edit')
    vals = strtrim(get(handles.v3, 'String'));
elseif strcmp(st,'popupmenu')
    contents = get(handles.v3,'String');
    vals = strtrim(char(contents(get(handles.v3, 'Value'))));
end

Keyword1 = handles.userdata.keyword.selected;
varName = sprintf('handles.userdata.%s.vars{3}', Keyword1);
Keyword2 = eval(varName);
vars = sprintf('handles.userdata.%s.%s', Keyword1, Keyword2);

type = eval(sprintf('handles.userdata.%s.args(3)', Keyword1));

if type <10
    strC = sprintf('%s=str2num(vals);',vars);
    eval(strC);
else
    strC = sprintf('%s=vals;',vars);
    eval(strC);
end

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes during object creation, after setting all properties.
function v3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

return;


function v4_Callback(hObject, eventdata, handles)
% hObject    handle to v4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of v4 as text
%        str2double(get(hObject,'String')) returns contents of v4 as a double
st = get(handles.v4, 'Style');

if strcmp(st,'edit')
    vals = strtrim(get(handles.v4, 'String'));
elseif strcmp(st,'popupmenu')
    contents = get(handles.v4,'String');
    vals = strtrim(char(contents(get(handles.v4, 'Value'))));
end

Keyword1 = handles.userdata.keyword.selected;
varName = sprintf('handles.userdata.%s.vars{4}', Keyword1);
Keyword2 = eval(varName);
vars = sprintf('handles.userdata.%s.%s', Keyword1, Keyword2);

type = eval(sprintf('handles.userdata.%s.args(4)', Keyword1));

if type <10
    strC = sprintf('%s=str2num(vals);',vars);
    eval(strC);
else
    strC = sprintf('%s=vals;',vars);
    eval(strC);
end

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes during object creation, after setting all properties.
function v4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


return;


function v5_Callback(hObject, eventdata, handles)
% hObject    handle to v5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of v5 as text
%        str2double(get(hObject,'String')) returns contents of v5 as a double
st = get(handles.v5, 'Style');

if strcmp(st,'edit')
    vals = strtrim(get(handles.v5, 'String'));
elseif strcmp(st,'popupmenu')
    contents = get(handles.v5,'String');
    vals = strtrim(char(contents(get(handles.v5, 'Value'))));
end

Keyword1 = handles.userdata.keyword.selected;
varName = sprintf('handles.userdata.%s.vars{5}', Keyword1);
Keyword2 = eval(varName);
vars = sprintf('handles.userdata.%s.%s', Keyword1, Keyword2);

type = eval(sprintf('handles.userdata.%s.args(5)', Keyword1));

if type <10
    strC = sprintf('%s=str2num(vals);',vars);
    eval(strC);
else
    strC = sprintf('%s=vals;',vars);
    eval(strC);
end

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes during object creation, after setting all properties.
function v5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


return;


function v6_Callback(hObject, eventdata, handles)
% hObject    handle to v6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of v6 as text
%        str2double(get(hObject,'String')) returns contents of v6 as a double
st = get(handles.v6, 'Style');

if strcmp(st,'edit')
    vals = strtrim(get(handles.v6, 'String'));
elseif strcmp(st,'popupmenu')
    contents = get(handles.v6,'String');
    vals = strtrim(char(contents(get(handles.v6, 'Value'))));
end

Keyword1 = handles.userdata.keyword.selected;
varName = sprintf('handles.userdata.%s.vars{6}', Keyword1);
Keyword2 = eval(varName);
vars = sprintf('handles.userdata.%s.%s', Keyword1, Keyword2);

type = eval(sprintf('handles.userdata.%s.args(6)', Keyword1));

if type <10
    strC = sprintf('%s=str2num(vals);',vars);
    eval(strC);
else
    strC = sprintf('%s=vals;',vars);
    eval(strC);
end

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes during object creation, after setting all properties.
function v6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


return;


function v7_Callback(hObject, eventdata, handles)
% hObject    handle to v7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of v7 as text
%        str2double(get(hObject,'String')) returns contents of v7 as a double
st = get(handles.v7, 'Style');

if strcmp(st,'edit')
    vals = strtrim(get(handles.v7, 'String'));
elseif strcmp(st,'popupmenu')
    contents = get(handles.v7,'String');
    vals = strtrim(char(contents(get(handles.v7, 'Value'))));
end

Keyword1 = handles.userdata.keyword.selected;
varName = sprintf('handles.userdata.%s.vars{7}', Keyword1);
Keyword2 = eval(varName);
vars = sprintf('handles.userdata.%s.%s', Keyword1, Keyword2);

type = eval(sprintf('handles.userdata.%s.args(7)', Keyword1));

if type <10
    strC = sprintf('%s=str2num(vals);',vars);
    eval(strC);
else
    strC = sprintf('%s=vals;',vars);
    eval(strC);
end

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes during object creation, after setting all properties.
function v7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


return;


function v8_Callback(hObject, eventdata, handles)
% hObject    handle to v8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of v8 as text
%        str2double(get(hObject,'String')) returns contents of v8 as a double
st = get(handles.v8, 'Style');

if strcmp(st,'edit')
    vals = strtrim(get(handles.v8, 'String'));
elseif strcmp(st,'popupmenu')
    contents = get(handles.v8,'String');
    vals = strtrim(char(contents(get(handles.v8, 'Value'))));
end

Keyword1 = handles.userdata.keyword.selected;
varName = sprintf('handles.userdata.%s.vars{8}', Keyword1);
Keyword2 = eval(varName);
vars = sprintf('handles.userdata.%s.%s', Keyword1, Keyword2);

type = eval(sprintf('handles.userdata.%s.args(8)', Keyword1));

if type <10
    strC = sprintf('%s=str2num(vals);',vars);
    eval(strC);
else
    strC = sprintf('%s=vals;',vars);
    eval(strC);
end

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes during object creation, after setting all properties.
function v8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


return;


function v9_Callback(hObject, eventdata, handles)
% hObject    handle to v9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of v9 as text
%        str2double(get(hObject,'String')) returns contents of v9 as a double
st = get(handles.v9, 'Style');

if strcmp(st,'edit')
    vals = strtrim(get(handles.v9, 'String'));
elseif strcmp(st,'popupmenu')
    contents = get(handles.v9,'String');
    vals = strtrim(char(contents(get(handles.v9, 'Value'))));
end

Keyword1 = handles.userdata.keyword.selected;
varName = sprintf('handles.userdata.%s.vars{9}', Keyword1);
Keyword2 = eval(varName);
vars = sprintf('handles.userdata.%s.%s', Keyword1, Keyword2);

type = eval(sprintf('handles.userdata.%s.args(9)', Keyword1));

if type <10
    strC = sprintf('%s=str2num(vals);',vars);
    eval(strC);
else
    strC = sprintf('%s=vals;',vars);
    eval(strC);
end

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes during object creation, after setting all properties.
function v9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


return;


function v10_Callback(hObject, eventdata, handles)
% hObject    handle to v10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of v10 as text
%        str2double(get(hObject,'String')) returns contents of v10 as a double
st = get(handles.v10, 'Style');

if strcmp(st,'edit')
    vals = strtrim(get(handles.v10, 'String'));
elseif strcmp(st,'popupmenu')
    contents = get(handles.v10,'String');
    vals = strtrim(char(contents(get(handles.v10, 'Value'))));
end

Keyword1 = handles.userdata.keyword.selected;
varName = sprintf('handles.userdata.%s.vars{10}', Keyword1);
Keyword2 = eval(varName);
vars = sprintf('handles.userdata.%s.%s', Keyword1, Keyword2);

type = eval(sprintf('handles.userdata.%s.args(10)', Keyword1));

if type <10
    strC = sprintf('%s=str2num(vals);',vars);
    eval(strC);
else
    strC = sprintf('%s=vals;',vars);
    eval(strC);
end

% Update handles structure
guidata(hObject, handles);

return;


% --- Executes during object creation, after setting all properties.
function v10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to v10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


return;


function eraseOptions(handles)

for i = 1:10
    strVT = sprintf('handles.vt%d',i);
    strV = sprintf('handles.v%d', i);
    strVC = sprintf('handles.vc%d', i);
    set(eval(strVT), 'Visible', 'off');
    set(eval(strV), 'Style', 'edit', 'Visible', 'off', 'String', '');
    set(eval(strVC), 'Visible', 'off');
end

set(handles.IOT_tag, 'Visible', 'off', 'String', 'Select Input','Position',handles.userdata.IOT_tag_pos);
set(handles.IO_tag, 'Visible', 'off','Position',handles.userdata.IO_tag_pos);
set(handles.IOL_tag, 'Visible', 'off', 'String', '','Position',handles.userdata.IOL_tag_pos);

set(handles.IOT2_tag, 'Visible', 'off', 'String', 'Select Input','Position',handles.userdata.IOT_tag_pos);
set(handles.IO2_tag, 'Visible', 'off','Position',handles.userdata.IO_tag_pos);
set(handles.IOL2_tag, 'Visible', 'off', 'String', '','Position',handles.userdata.IOL_tag_pos);

return;


% --- Executes on button press in vc1.
function vc1_Callback(hObject, eventdata, handles)
% hObject    handle to vc1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = sprintf('handles.userdata.%s.args(1)', handles.userdata.keyword.selected);
vals = eval(str);

if vals == 200
    pathname = uigetdir(pwd, 'Select a directory');
    if pathname ~= 0
        set(handles.v1, 'String', pathname);
    end
elseif vals == 100
    [filename, pathname] = uigetfile({'*.*','All Files (*.*)'}, 'Pick a file');
    if filename ~= 0
        set(handles.v1, 'String', fullfile(pathname, filename));
    end
end

v1_Callback(hObject, eventdata, handles);

return;

% --- Executes during object creation, after setting all properties.
function vc1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vc1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

return;

% --- Executes on button press in vc2.
function vc2_Callback(hObject, eventdata, handles)
% hObject    handle to vc2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = sprintf('handles.userdata.%s.args(2)', handles.userdata.keyword.selected);
vals = eval(str);

if vals == 200
    pathname = uigetdir(pwd, 'Select a directory');
    if pathname ~= 0
        set(handles.v2, 'String', pathname);
    end
elseif vals == 100
    [filename, pathname] = uigetfile({'*.*','All Files (*.*)'}, 'Pick a file');
    if filename ~= 0
        set(handles.v2, 'String', fullfile(pathname, filename));
    end
end

v2_Callback(hObject, eventdata, handles);

return;

% --- Executes on button press in vc3.
function vc3_Callback(hObject, eventdata, handles)
% hObject    handle to vc3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = sprintf('handles.userdata.%s.args(3)', handles.userdata.keyword.selected);
vals = eval(str);

if vals == 200
    pathname = uigetdir(pwd, 'Select a directory');
    if pathname ~= 0
        set(handles.v3, 'String', pathname);
    end
elseif vals == 100
    [filename, pathname] = uigetfile({'*.*','All Files (*.*)'}, 'Pick a file');
    if filename ~= 0
        set(handles.v3, 'String', fullfile(pathname, filename));
    end
end

v3_Callback(hObject, eventdata, handles);

return;

% --- Executes on button press in vc4.
function vc4_Callback(hObject, eventdata, handles)
% hObject    handle to vc4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = sprintf('handles.userdata.%s.args(4)', handles.userdata.keyword.selected);
vals = eval(str);

if vals == 200
    pathname = uigetdir(pwd, 'Select a directory');
    if pathname ~= 0
        set(handles.v4, 'String', pathname);
    end
elseif vals == 100
    [filename, pathname] = uigetfile({'*.*','All Files (*.*)'}, 'Pick a file');
    if filename ~= 0
        set(handles.v4, 'String', fullfile(pathname, filename));
    end
end

v4_Callback(hObject, eventdata, handles);

return;

% --- Executes on button press in vc5.
function vc5_Callback(hObject, eventdata, handles)
% hObject    handle to vc5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = sprintf('handles.userdata.%s.args(5)', handles.userdata.keyword.selected);
vals = eval(str);

if vals == 200
    pathname = uigetdir(pwd, 'Select a directory');
    if pathname ~= 0
        set(handles.v5, 'String', pathname);
    end
elseif vals == 100
    [filename, pathname] = uigetfile({'*.*','All Files (*.*)'}, 'Pick a file');
    if filename ~= 0
        set(handles.v5, 'String', fullfile(pathname, filename));
    end
end

v5_Callback(hObject, eventdata, handles);

return;

% --- Executes on button press in vc6.
function vc6_Callback(hObject, eventdata, handles)
% hObject    handle to vc6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = sprintf('handles.userdata.%s.args(6)', handles.userdata.keyword.selected);
vals = eval(str);

if vals == 200
    pathname = uigetdir(pwd, 'Select a directory');
    if pathname ~= 0
        set(handles.v6, 'String', pathname);
    end
elseif vals == 100
    [filename, pathname] = uigetfile({'*.*','All Files (*.*)'}, 'Pick a file');
    if filename ~= 0
        set(handles.v6, 'String', fullfile(pathname, filename));
    end
end

v6_Callback(hObject, eventdata, handles);

return;

% --- Executes on button press in vc7.
function vc7_Callback(hObject, eventdata, handles)
% hObject    handle to vc7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = sprintf('handles.userdata.%s.args(7)', handles.userdata.keyword.selected);
vals = eval(str);

if vals == 200
    pathname = uigetdir(pwd, 'Select a directory');
    if pathname ~= 0
        set(handles.v7, 'String', pathname);
    end
elseif vals == 100
    [filename, pathname] = uigetfile({'*.*','All Files (*.*)'}, 'Pick a file');
    if filename ~= 0
        set(handles.v7, 'String', fullfile(pathname, filename));
    end
end

v7_Callback(hObject, eventdata, handles);

return;

% --- Executes on button press in vc8.
function vc8_Callback(hObject, eventdata, handles)
% hObject    handle to vc8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = sprintf('handles.userdata.%s.args(8)', handles.userdata.keyword.selected);
vals = eval(str);

if vals == 200
    pathname = uigetdir(pwd, 'Select a directory');
    if pathname ~= 0
        set(handles.v8, 'String', pathname);
    end
elseif vals == 100
    [filename, pathname] = uigetfile({'*.*','All Files (*.*)'}, 'Pick a file');
    if filename ~= 0
        set(handles.v8, 'String', fullfile(pathname, filename));
    end
end

v8_Callback(hObject, eventdata, handles);

return;

% --- Executes on button press in vc9.
function vc9_Callback(hObject, eventdata, handles)
% hObject    handle to vc9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = sprintf('handles.userdata.%s.args(9)', handles.userdata.keyword.selected);
vals = eval(str);

if vals == 200
    pathname = uigetdir(pwd, 'Select a directory');
    if pathname ~= 0
        set(handles.v9, 'String', pathname);
    end
elseif vals == 100
    [filename, pathname] = uigetfile({'*.*','All Files (*.*)'}, 'Pick a file');
    if filename ~= 0
        set(handles.v9, 'String', fullfile(pathname, filename));
    end
end

v9_Callback(hObject, eventdata, handles);

return;

% --- Executes on button press in vc10.
function vc10_Callback(hObject, eventdata, handles)
% hObject    handle to vc10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

str = sprintf('handles.userdata.%s.args(10)', handles.userdata.keyword.selected);
vals = eval(str);

if vals == 200
    pathname = uigetdir(pwd, 'Select a directory');
    if pathname ~= 0
        set(handles.v10, 'String', pathname);
    end
elseif vals == 100
    [filename, pathname] = uigetfile({'*.*','All Files (*.*)'}, 'Pick a file');
    if filename ~= 0
        set(handles.v10, 'String', fullfile(pathname, filename));
    end
end

v10_Callback(hObject, eventdata, handles);

return;


% --- Executes during object creation, after setting all properties.
function LoadTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LoadTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

return;


% --- Executes during object creation, after setting all properties.
function SaveTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SaveTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

return;

% --- Executes during object creation, after setting all properties.
function RunTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RunTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

return;

% --- Executes during object creation, after setting all properties.
function CancelTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to CancelTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

return;

% --- Executes during object creation, after setting all properties.
function IO_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IO_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

return;


% --- Executes during object creation, after setting all properties.
function IOT_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IOT_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

return;


% --- Executes on button press in IO2_tag.
function IO2_tag_Callback(hObject, eventdata, handles)
% hObject    handle to IO2_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

switch handles.userdata.keyword.selected
    case 'ParamCALD'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.ParamCALD.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'ParamPDM'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.ParamPDM.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'ParaQuad'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.ParamQuad.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'ParamConf'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.ParamConf.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'ExpLSF'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.ExpLSF.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'ExpIRF'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.ExpIRF.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'AligSHREC'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.AligSHREC.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'AligFOE'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.AligFOE.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'AligLandmark'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.AligLandmark.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'ExpAlig'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.ExpAlig.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'PCA'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.PCA.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 't_map'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.t_map.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'Import'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.Import.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'res_t_map'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.res_t_map.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'off');
    case 'res_PCA'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.res_PCA.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'off');
    case 'bim2gipl'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.bim2gipl.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'fix2gipl'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.fix2gipl.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'gipl2bim'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.gipl2bim.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'smo2surf_para_meta'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.smo2surf_para_meta.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'surf_para_meta2smo'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.surf_para_meta2smo.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'des2meta_coef'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.des2meta_coef.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'meta_coef2des'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.meta_coef2des.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'ellalign_meta_coef2reg'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.ellalign_meta_coef2reg.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'reg2procalign_meta'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.reg2procalign_meta.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'InHouse_Fix'
%        setappdata(0, 'UseNativeSystemDialogs', 1);
        [filename, pathname, filterindex] = uigetfile(handles.userdata.InHouse_Fix.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'PDM_Fix'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.PDM_Fix.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'DisplayObjs'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.DisplayObjs.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'AverageObjs'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.AverageObjs.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    case 'ScaleObjs'
        [filename, pathname, filterindex] = uigetfile(handles.userdata.ScaleObjs.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
    otherwise
        [filename, pathname, filterindex] = uigetfile(handles.userdata.others.inFilter, ...
        'Pick Input Objects', 'MultiSelect', 'on');
end

if filterindex ~= 0
    if (class(filename) == 'cell') 
        d = size(filename,2);
        for i = 1:d
            files{i} = [pathname filename{i}];
        end
        set(handles.IOL2_tag, 'String', files);

    else
        d = size(filename,1);
        for i = 1:d
            files{i} = [pathname filename(i,:)];
        end
        set(handles.IOL2_tag, 'String', files);        
    end
    handles.userdata.inObjs2 = files;
end
% Update handles structure
guidata(hObject, handles);

return;


function IOL2_tag_Callback(hObject, eventdata, handles)
% hObject    handle to IOL2_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of IOL2_tag as text
%        str2double(get(hObject,'String')) returns contents of IOL2_tag as a double
handles.userdata.inObjs2 = get(hObject, 'String');

% Update handles structure
guidata(hObject, handles);
return;


% --- Executes during object creation, after setting all properties.
function IOL2_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IOL2_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


return;


function updateVariables(hObject, eventdata, handles)

return;
