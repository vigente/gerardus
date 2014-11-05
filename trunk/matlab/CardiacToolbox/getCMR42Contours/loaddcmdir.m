function [dcmSeries, dcmPatient, SeriesList, SeriesListNumInfo, ImageInfos, SeriesInfos] = loaddcmdir(dicomdirpath)
%function [dcmSeries, dcmPatient, SeriesList, SeriesListNumInfo, ImageInfos, SeriesInfos] = loaddcmdir(dicomdirpath)
%
% LOADDICOMDIR reads metadata from DICOMDIR file and shows list of
% the different DICOM SERIES in Pop-up menus. After selecting a dicom 
% series the program creates the dcmSeries structure containig the 
% related file names and path.
% 
% Input: (OPTIONAL)
%   dicomdirpath  - the path where the DICOMDIR file is located                       
%   
% Output:
%   dcmSeries.Images        - cell array containing the file names of the images
%   dcmSeries.Path          - directory of the image files
%   dcmSeries.dicomdirPath  - directory of DICOMDIR file
%
%   dcmPatient              - the full dicom structure from DICOMDIR
%   SeriesList              - the full list of the available series from DICOMDIR
%                               (cell array)
%   SeriesListNumInfo       - structure containing the patients, study and series 
%                               number relating to the series list numbers
%   ImageInfos              - For all DirectoryRecordSequence field of the DICOMDIR structure ...
%                           the following data are listed  ReferencedFileID, InstanceNumber, ...
%                           ImageType, ReferencedSOPInstanceUIDInFile, ...
%                           SeriesInstanceUID, CurrentSeries, SeriesNumber
%                           (structure)
%   SeriesInfos             - SeriesInstanceUID and SeriesNumber data for
%                           all series (structure)
%
% Hints: 
% 
% Reading DICOMDIR folder using the directory browser GUI:
% dcmSeries = loaddcmdir;
% 
% Generating 3D array by selecting a dicom series:
% [imaVOL, scaninfo, dcminfo] = loaddcm(loaddcmdir);
% 
% To preview the 3D array you can use the orthogonalslicer tool: (http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=5934&objectType=file)
% 
% orthogonalslicer(imaVOL,scaninfo.pixsize,'gray');
% 
% To check the dicom header:
% openvar('dcminfo'); 
%
%
% Matlab library function for MIA_gui utility. 
% University of Debrecen, PET Center/LB 2008

%try
    hm = [];
    dcmSeries = []; dcmPatient = []; SeriesList = [];  SeriesListNumInfo = []; ImageInfos = []; SeriesInfos = [];
    ImgPrcTlbx = ver('images');
    ImgPrcTlbxVersion = str2double(ImgPrcTlbx.Version(1));
    if ImgPrcTlbxVersion < 4
        disp('The required minimum version of ImgProc. Toolbox is 4!');
        return; 
    end
    
	if nargin == 0 | strcmp(dicomdirpath,'')  
         dicomdirpath = uigetdir('','Select DICOMDIR folder');
         if dicomdirpath == 0
             return;
         end
    end
    filename1mat = fullfile(dicomdirpath,'dicomdir.mat');  
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
%{
        % NOTE commented out by Tasos
        disp(' ');                     
        disp('load dicomdir.mat');
        disp(' ');
%}        
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
            CurrentItem = getfield(DirRecords,DirRecordsFieldNames{i});
            CurrentItemType = getfield(CurrentItem,'DirectoryRecordType');
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
                    SeriesList{SeriesListNum} = [Pname,',  ',Modality,',  SeriesDescr: ', SeriesDescription, ...
                        ',  StudyDesc: ', StudyDescription, ',  Number of Images: ',NumOfImages];
                    SeriesListNumInfo(SeriesListNum).PatientNum = CurrentPatient;
                    SeriesListNumInfo(SeriesListNum).StudyNum = CurrentStudy;
                    SeriesListNumInfo(SeriesListNum).SeriesNum = CurrentSeries;
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
                if isfield(CurrentItem,'SeriesInstanceUID')
                    SeriesInfos{CurrentSeries,1} = CurrentItem.SeriesInstanceUID;
                else
                    SeriesInfos{CurrentSeries,1} = '';
                end
                if isfield(CurrentItem,'SeriesNumber')
                    SeriesInfos{CurrentSeries,2} = CurrentItem.SeriesNumber;
                else
                    SeriesInfos{CurrentSeries,2} = 0;
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
                ImageInfos{i,1} = [pathstr,filesep,dcmfname];
                if isfield(CurrentItem,'InstanceNumber')
                    ImageInfos{i,2} = CurrentItem.InstanceNumber;
                else
                    ImageInfos{i,2} = 0;
                end
                if isfield(CurrentItem,'ImageType')
                    ImageInfos{i,3} = CurrentItem.ImageType;
                else
                    ImageInfos{i,3} = '';
                end
                ImageInfos{i,4} = CurrentItem.ReferencedSOPInstanceUIDInFile;
                ImageInfos{i,5} = SeriesInfos{CurrentSeries,1};
                ImageInfos{i,6} = CurrentSeries;
                ImageInfos{i,7} = SeriesInfos{CurrentSeries,2};
                %%%%
            end
            waitbar(i/NumOfField,waitbar_h);
            drawnow;
        end
        close(waitbar_h);
        % create summary about the last read series for popupmenu
        SeriesListNum = SeriesListNum + 1;
        Pname = dcmPatient(CurrentPatient,1).PatientName;
        Modality = dcmPatient(CurrentPatient,1).Study(CurrentStudy,1).Series(CurrentSeries,1).Modality;
        SeriesDescription = dcmPatient(CurrentPatient,1).Study(CurrentStudy,1). ...
            Series(CurrentSeries,1).SeriesDescription;
        NumOfImages = num2str(CurrentImage);
        SeriesList{SeriesListNum} = [Pname,',  ',Modality,',  SeriesDesc: ', SeriesDescription, ...
            ',  StudyDesc: ', StudyDescription,',  Number of Images: ',NumOfImages];
        SeriesListNumInfo(SeriesListNum).PatientNum = CurrentPatient;
        SeriesListNumInfo(SeriesListNum).StudyNum = CurrentStudy;
        SeriesListNumInfo(SeriesListNum).SeriesNum = CurrentSeries;
        
        try
            dcmdir.SeriesListNumInfo = SeriesListNumInfo;
            dcmdir.SeriesList = SeriesList;
            dcmdir.dcmPatient = dcmPatient;
            dcmdir.dcmSeries = dcmSeries;
            save(filename1mat,'dcmdir');
        catch
            % Unable to write file 'filename1mat' permission denied.
        end
        
    end
    
    % create popupmenu for dicomseries
    global dcmdirlistVal;
    dcmdirlistVal = 0;
% NOTE +500 in next two lines added by tasos to enlarge window
    dcmdirlistfh = figure('menubar','none','NumberTitle','off','name','DICOMDIR info','position',[250   400   760+500   520+500]);
    lbh = uicontrol('Style','listbox','Position',[10 60 720+500 440+500],'tag','dcmdirlist_popupmenu');
    set(lbh,'string',SeriesList);
    OKCallback = 'global dcmdirlistVal; dcmdirlistVal = get(findobj(''tag'',''dcmdirlist_popupmenu''),''value''); delete(findobj(''name'',''DICOMDIR info''));';
    CancelCallback = 'delete(findobj(''name'',''DICOMDIR info''));';
    OK_h = uicontrol('Style', 'pushbutton', 'String', 'OK','Position', [440 10 80 30], 'Callback', OKCallback);
    Cancel_h = uicontrol('Style', 'pushbutton', 'String', 'Cancel','Position', [340 10 80 30], 'Callback', CancelCallback);
    
    uiwait(dcmdirlistfh);
    
    if dcmdirlistVal == 0;
        disp('No DICOM Series was selected!');
        return;
    end
    
    % create the outputs
    dcmSeriesPath = dcmPatient(SeriesListNumInfo(dcmdirlistVal).PatientNum). ...
        Study(SeriesListNumInfo(dcmdirlistVal).StudyNum). ...
        Series(SeriesListNumInfo(dcmdirlistVal).SeriesNum). ...
        ImagePath;
    dcmSeries.Path = [dicomdirpath,filesep];
    dcmSeries.dicomdirpath = [dicomdirpath, filesep, dcmSeriesPath, filesep];
    %dcmSeries.Path = [dicomdirpath, filesep, dcmSeriesPath, filesep];
    %dcmSeries.dicomdirpath = dicomdirpath;
    dcmSeries.Images = dcmPatient(SeriesListNumInfo(dcmdirlistVal).PatientNum). ...
        Study(SeriesListNumInfo(dcmdirlistVal).StudyNum). ...
        Series(SeriesListNumInfo(dcmdirlistVal).SeriesNum). ...
        ImageNames;
    
% catch %in case of any error
%     ErrorOnDicomDirOpening = lasterr
%     if isobject(hm)
%         delete(hm);
%     end
% end