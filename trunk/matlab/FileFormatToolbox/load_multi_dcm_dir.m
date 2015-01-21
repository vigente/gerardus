function [seriesData, dicomdirpath] = load_multi_dcm_dir(dicomdirpath)
% LOAD_MULTI_DCM_DIR  Read metadata from DICOMDIR file.
%
% LOAD_MULTI_DCM_DIR reads metadata from DICOMDIR file and shows list of
% the different DICOM SERIES in Pop-up menus. After selecting a dicom 
% series the program creates the dcmSeries structure containing the 
% related file names and path.
%
% DICOMDIR is a standard binary file provided with DICOM data, and has
% information such as paths to files, and patient, study, series IDs, etc.
% 
% SERIESDATA = LOAD_MULTI_DCM_DIR(DICOMDIRPATH)
%
%   DICOMDIRPATH is a string with the path to the directory containing the
%   DICOMDIR file. If not provided, it opens a file browser to select the
%   directory graphically.
%
%   SERIESDATA is a structure containing image information for series 
%   selected by the user in the GUI.

% Original Author: Laszlo Balkay, University of Debrecen, PET Center/LB 2008
%      Available from: http://uk.mathworks.com/matlabcentral/fileexchange/7926-dicomdir-reader
%      under BSD 2-Clause License.
%
% Modified by Christopher Kelly <christopher.kelly28@googlemail.com> and
% Benjamin Villard <b.016434@gmail.com>, University of Oxford.
% Copyright © 2014-2015 University of Oxford
% Version: 0.2.0
% $Rev$
% $Date$
%
% University of Oxford means the Chancellor, Masters and Scholars of
% the University of Oxford, having an administrative office at
% Wellington Square, Oxford OX1 2JD, UK. 
%
% 10 Dec 2014 - Interface changed to give preview of image. 
%             - Also allows the selection of multiple series. 
%             - Fewer output arguments, now the only output contains the
%               information of the series selected
%
% This file is distributed as a derivative work of a third-party function
% with project Gerardus.
%
% http://code.google.com/p/gerardus/

% check arguments
narginchk(0, 1);
nargoutchk(0, 2);

hm = [];
dcmPatient = [];
ImgPrcTlbx = ver('images');
ImgPrcTlbxVersion = str2double(ImgPrcTlbx.Version(1));
if ImgPrcTlbxVersion < 4
    disp('The required minimum version of ImgProc. Toolbox is 4!');
    return;
end

if ((nargin == 0) || strcmp(dicomdirpath,''))
    dicomdirpath = uigetdir('','Select DICOMDIR folder');
    if dicomdirpath == 0
        return;
    end
end
filename1mat = fullfile(dicomdirpath,'dicom.mat');
dirres1mat = dir(filename1mat);
filename1 = fullfile(dicomdirpath,'dicomdir');
dirres1 = dir(filename1);
filename2 = fullfile(dicomdirpath,'DICOMDIR');
dirres2 = dir(filename2);
if ~isempty(dirres1mat)
    % if the dicomdir mat structure was saved on the first time
    % open it. It is much faster than reading the DICOMDIR dicom file
    load(filename1mat);
    SeriesListNumInfo = dcmdir.SeriesListNumInfo;
    SeriesList = dcmdir.SeriesList;
    dcmPatient = dcmdir.dcmPatient;
    dcmSeries = dcmdir.dcmSeries;
    clear dcmdir;
    disp(' ');
    disp('load dicomdir.mat');
    disp(' ');
else
    if ~isempty(dirres1)
        dcmdir_path = filename1;
    elseif ~isempty(dirres2)
        dcmdir_path = filename2;
    else
        hm = msgbox('No DICOMDIR in the selected folder !','MIA Info','warn' );
        disp('No DICOMDIR in the selected folder !')
        return;
    end
    
    
    % Show up message window and set the cursor type to "watch"
    hm = msgbox('Reading the DICOMDIR structure. It takes time. Please wait!','MIA Info' );
    set(hm, 'tag', 'Msgbox_MIA Info');
    SetData=setptr('watch');set(hm,SetData{:});
    hmc = (get(hm,'children'));
    set(hmc(2),'enable','inactive');
    pause(1);
    
    % open DICOMDIR file & read DirectoryRecordSequence structure
    dcmhdr = dicominfo(dcmdir_path);
    delete(hm);
    
    % Read the DirectoryRecordSequence structure until end,
    % and read the all included DirectoryRecordType item sequentially
    DirRecords = dcmhdr.DirectoryRecordSequence;
    DirRecordsFieldNames = fieldnames(DirRecords);
    NumOfField = size(DirRecordsFieldNames,1);
    CurrentPatient = 0;
    SeriesListNum = 0;
    waitbar_h = waitbar(0,'Scanning the DICOM files....');
    for i = 1: NumOfField
        CurrentItem = DirRecords.(DirRecordsFieldNames{i});
        CurrentItemType = CurrentItem.DirectoryRecordType;
        if strcmp(CurrentItemType,'PATIENT')
            CurrentPatient = CurrentPatient + 1;
            CurrentStudy = 0;
            if ImgPrcTlbxVersion > 4
                if isfield(CurrentItem,'PatientName')
                    dcmPatient(CurrentPatient,1).PatientName = ...
                        CurrentItem.PatientName.FamilyName;
                else
                    dcmPatient(CurrentPatient,1).PatientName = 'No value available';
                end
            elseif ImgPrcTlbxVersion == 4
                if isfield(CurrentItem,'PatientsName')
                    dcmPatient(CurrentPatient,1).PatientName = ...
                        CurrentItem.PatientsName.FamilyName;
                else
                    dcmPatient(CurrentPatient,1).PatientName = 'No value available';
                end
            end
        elseif strcmp(CurrentItemType,'STUDY')
            CurrentStudy = CurrentStudy + 1;
            CurrentSeries = 0;
            if isfield(CurrentItem,'StudyDescription')
                dcmPatient(CurrentPatient,1).Study(CurrentStudy,1).StudyDescription = ...
                    CurrentItem.StudyDescription;
            else
                dcmPatient(CurrentPatient,1).Study(CurrentStudy,1).StudyDescription = '';
            end
        elseif strcmp(CurrentItemType,'SERIES')
            if CurrentSeries > 0 % create summary about the previous read series for popupmenu
                SeriesListNum = SeriesListNum + 1;
                Pname = dcmPatient(CurrentPatient,1).PatientName;
                Modality = dcmPatient(CurrentPatient,1).Study(CurrentStudy,1).Series(CurrentSeries,1).Modality;
                SeriesDescription = dcmPatient(CurrentPatient,1).Study(CurrentStudy,1). ...
                    Series(CurrentSeries,1).SeriesDescription;
                StudyDescription = dcmPatient(CurrentPatient,1).Study(CurrentStudy,1).StudyDescription;
                NumOfImages = num2str(CurrentImage);
            end
            CurrentSeries = CurrentSeries + 1;
            CurrentImage = 0;
            dcmPatient(CurrentPatient,1).Study(CurrentStudy,1).Series(CurrentSeries,1).Modality = ...
                CurrentItem.Modality;
            if isfield(CurrentItem,'SeriesDescription')
                dcmPatient(CurrentPatient,1).Study(CurrentStudy,1).Series(CurrentSeries,1). ...
                    SeriesDescription = CurrentItem.SeriesDescription;
            else
                dcmPatient(CurrentPatient,1).Study(CurrentStudy,1).Series(CurrentSeries,1). ...
                    SeriesDescription = '';
            end
            %%%
        elseif strcmp(CurrentItemType,'IMAGE')
            CurrentImage = CurrentImage + 1;
            [pathstr, dcmfname, ext] = fileparts(CurrentItem.ReferencedFileID);
            if CurrentImage == 1
                dcmPatient(CurrentPatient,1).Study(CurrentStudy,1).Series(CurrentSeries,1). ...
                    ImagePath = pathstr;
            end
            if ispc
                dcmPatient(CurrentPatient,1).Study(CurrentStudy,1).Series(CurrentSeries,1). ...
                    ImageNames{CurrentImage,1} = CurrentItem.ReferencedFileID;
            else
                dcmPatient(CurrentPatient,1).Study(CurrentStudy,1).Series(CurrentSeries,1). ...
                    ImageNames{CurrentImage,1} = strrep(CurrentItem.ReferencedFileID,'\','/');
            end
            %ImageNames{CurrentImage,1} = dcmfname;
            dcmPatient(CurrentPatient,1).Study(CurrentStudy,1).Series(CurrentSeries,1). ...
                ReferencedSOPInstanceUIDInFile{CurrentImage,1} = CurrentItem.ReferencedSOPInstanceUIDInFile;
            %%%
        end
        waitbar(i/NumOfField,waitbar_h);
        drawnow;
    end
    close(waitbar_h);
    
    try
        dcmdir.dcmPatient = dcmPatient;
        save(filename1mat,dcmdir);
    catch
        % Unable to write file 'filename1mat' permission denied.
    end
    
end

global dcmVal;
dcmVal = 0;

if size(dcmPatient,1) > 1
    selectPatientFigure = figure('Units','normalized','menubar','none',...
        'NumberTitle','off','name','selectPatient',...
        'position',[2/8, 3/8, 4/8, 2/8]);
    uicontrol('Parent',selectPatientFigure,'Units','normalized','Style','text',...
        'String','Select which patient you are analysing?','Position',[0.1 0.7 0.8 0.2]);
    String  = sprintf('%s|', dcmPatient.('PatientName'));
    String(end) = [];
    lbh = uicontrol('Parent',selectPatientFigure,'Units','normalized','Style','popupmenu',...
        'String',String,'Position',[0.1 0.4 0.8 0.2],'tag','patientList');
    okCloseCallback = 'global dcmVal; dcmVal = get(findobj(''tag'',''patientList''),''value''); delete(findobj(''name'',''selectPatient''));';
    uicontrol('Parent',selectPatientFigure,'Units','normalized','Style','pushbutton','String','OK','Position',[0.1 0.1 0.8 0.2],'Callback',okCloseCallback);
    uiwait(selectPatientFigure);
else
    dcmVal = 1;
end

dcmPatient = dcmPatient(dcmVal,:);

clearvars -except dcmPatient dicomdirpath

set(0,'DefaultFigureWindowStyle','normal');

screensize = get(0,'ScreenSize');
dcmdirlistfh = figure('Units','normalized','menubar','none','NumberTitle','off',...
    'name','DICOMDIR info','position',[0.1, 0.1, 0.8, 0.8]);

lbh = uicontrol('Parent',dcmdirlistfh,'Units','normalized','Style','listbox',...
    'Position',[0.05, 0.09, 0.45, 0.88 ],'tag','dcmdirlist_popupmenu');

[seriesDescriptions,order] = sort({dcmPatient.Study.Series.SeriesDescription});
dcmPatient.Study.Series = dcmPatient.Study.Series(order);

set(lbh,'string',seriesDescriptions,'max',1.1);

global dcmdirlistVal;
dcmdirlistVal = 0;

OKCallback = 'global a; a = 1; global dcmdirlistVal; dcmdirlistVal = get(findobj(''tag'',''dcmdirlist_popupmenu''),''value''); delete(findobj(''name'',''DICOMDIR info'')); return;';
CancelCallback = 'global a; a = 1; delete(findobj(''name'',''DICOMDIR info''));';
OK_h = uicontrol('Parent',dcmdirlistfh,'Units','normalized','Style','pushbutton',...
    'String','OK','Position',[0.15, 0.03, 0.2, 0.03],'Callback',OKCallback);
Cancel_h = uicontrol('Parent',dcmdirlistfh,'Units','normalized','Style','pushbutton',...
    'String','Cancel','Position',[0.65, 0.03, 0.2, 0.03],'Callback',CancelCallback);

global a
a = 0;
axes1 = axes('Parent',dcmdirlistfh,'Units','normalized','Position',[0.55, 0.55, 0.4, 0.4]); axis square;
while a == 0
    dcmdirlistVal = get(findobj('tag','dcmdirlist_popupmenu'),'value');
    pause(0.1);
    if length(dcmdirlistVal) >= 1.1
        continue
    end
    try
        imageCurrent = double(uint16(dicomread(fullfile(dicomdirpath,dcmPatient.Study.Series(dcmdirlistVal).ImageNames{1}))));
        imagesc(imageCurrent); axis square; colormap gray; axis off;
        pause(0.1);
        
        strDetails = sprintf('File Location : %s \n\n Series Description : %s \n\n Number of Frames : %s \n\n ',...
            num2str(fullfile(dicomdirpath,dcmPatient.Study.Series(dcmdirlistVal).ImageNames{1})),...
            num2str(dcmPatient.Study.Series(dcmdirlistVal).SeriesDescription),...
            num2str(numel(dcmPatient.Study.Series(dcmdirlistVal).ImageNames)));
        
        uicontrol('Parent',dcmdirlistfh,'Units','normalized','Style','text',...
            'Position',[0.55, 0.09, 0.4, 0.4],'String',strDetails,'FontWeight','bold');
    end
end

close all;

if dcmdirlistVal == 0;
    disp('No DICOM Series was selected!');
    return;
end

seriesData = dcmPatient.Study.Series(dcmdirlistVal);

end
