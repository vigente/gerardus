function [imaVOL, scaninfo, dcminfo] = loaddcm(FileNames,dir_path)
%function [imaVOL, scaninfo, dcminfo] = loaddcm(FileNames,dir_path)
%
% Inputs: (OPTIONAL)
%   FileNames - cell array containing the file names of images
%   dir_path  - directory of the image files
%
%   Without inputs the user can multiselect the desired files
%   using the uigetfiles GUI  
%   
% Outputs:
%   imaVOL    - 3D int16 array of the images
%   scaninfo  - short structure including the most important 
%               info for the volume
%   dcminfo   - dicom header structure relating to
%               the last opened images file /FileNames(end)/ 
%
% Matlab library function for MIA_gui utility. 
% University of Debrecen, PET Center/LB 2003



if isempty(whos('global','mia_dicom_scaled'))
% Ez a loaddcm fv-nek �ll�tja be, hogy az eredeti pixel�rt�kek vagy a SUV sk�l�zott
% �rt�kek olvas�djanak be
        mia_dicom_scaled = 1; % ha nincs ilyen global v�ltoz�, akkor SUV sk�la lesz sz�molva
else
    global mia_dicom_scaled;
end

try
    hm = []; p = [];
	if nargin == 0
         [FilesSelected, dir_path] = uigetfiles('*.dcm','Select DICOM file');
          if isempty(FilesSelected);
              imaVOL = [];scaninfo = [];
              return;
          end
          FileNames = sortrows(FilesSelected');
          filename = [dir_path,char(FileNames(1))];
          num_of_files = size(FileNames,1);
		  for i=1:num_of_files
            filelist(i).name = char(FileNames(i));
		  end
		  filelist = filelist';
	elseif nargin == 1
          FilesSelected = FileNames.Images;
          dir_path = FileNames.Path;
          FileNames = sortrows(FilesSelected);
          filename = [dir_path,char(FileNames(1))];
          num_of_files = size(FileNames,1);
		  for i=1:num_of_files
            filelist(i).name = char(FileNames(i));
		  end
		  filelist = filelist';
	elseif nargin == 2
          num_of_files = size(FileNames,1);
          filename = [dir_path,char(FileNames(1))];
	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% open DICOM file & get the needed image parameters (scaninfo) 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	hm = msgbox('Loading the Dicom volume...','Dicom Info' );
	SetData=setptr('watch');set(hm,SetData{:});
	hmc = (get(hm,'children'));
	set(hmc(2),'enable','inactive');
    pause(1);
	
    PhilipsDcmDyn = 0;
	dcminfo = dicominfo(filename);
	%
	% get info from dcm file and save them in the scaninfo structure 
	%
    ImgPrcTlbx = ver('images');
    ImgPrcTlbxVersion = str2double(ImgPrcTlbx.Version);
    if ImgPrcTlbxVersion < 4.2
		if isfield(dcminfo.PatientsName,'GivenName')
            scaninfo.pnm	    = [dcminfo.PatientsName.GivenName,' ',dcminfo.PatientsName.FamilyName];	
		else
            scaninfo.pnm	    = dcminfo.PatientsName.FamilyName;	
		end
    elseif ImgPrcTlbxVersion >= 4.2
        if isfield(dcminfo.PatientName,'GivenName')
            scaninfo.pnm	    = [dcminfo.PatientName.GivenName,' ',dcminfo.PatientName.FamilyName];	
		else
            scaninfo.pnm	    = dcminfo.PatientName.FamilyName;	
		end
    end
	scaninfo.brn        = [];
	if isfield(dcminfo,'PatientID')
        scaninfo.rid	    = dcminfo.PatientID;	
	else
        scaninfo.rid	    = [];
	end
	scaninfo.rin	    = [];
	if isfield(dcminfo,'AcquisitionDate')
		if ~isempty(dcminfo.AcquisitionDate)
			scaninfo.daty	    = dcminfo.AcquisitionDate(1:4);
			scaninfo.datm	    = dcminfo.AcquisitionDate(5:6);
			scaninfo.datd       = dcminfo.AcquisitionDate(7:8);
		end
	elseif isfield(dcminfo,'StudyDate')
        if ~isempty(dcminfo.StudyDate)
            scaninfo.daty	    = dcminfo.StudyDate(1:4);
			scaninfo.datm	    = dcminfo.StudyDate(5:6);
			scaninfo.datd       = dcminfo.StudyDate(7:8);
        end
	else
        scaninfo.daty	    = [];
		scaninfo.datm	    = [];
		scaninfo.datd       = [];
	end
	if isfield(dcminfo,'AcquisitionTime')
		if ~isempty(dcminfo.AcquisitionTime)
			scaninfo.timh	    = dcminfo.AcquisitionTime(1:2);
			scaninfo.timm	    = dcminfo.AcquisitionTime(3:4);
			scaninfo.tims	    = dcminfo.AcquisitionTime(5:6);
		end
	elseif isfield(dcminfo,'StudyTime')
        if ~isempty(dcminfo.StudyTime)
            scaninfo.timh	    = dcminfo.StudyTime(1:2);
			scaninfo.timm	    = dcminfo.StudyTime(3:4);
			scaninfo.tims	    = dcminfo.StudyTime(5:6);
        end
	else
        scaninfo.timh	    = [];
		scaninfo.timm	    = [];
		scaninfo.tims	    = [];
	end
	scaninfo.mtm        = [];
	if isfield(dcminfo,'RadiopharmaceuticalInformationSequence')
        if isfield(dcminfo.RadiopharmaceuticalInformationSequence.Item_1,'Radiopharmaceutical')
            scaninfo.iso 	= dcminfo.RadiopharmaceuticalInformationSequence.Item_1.Radiopharmaceutical;
        else
            scaninfo.iso 	= [];
        end
	else
        scaninfo.iso 	= [];
	end
	scaninfo.half       = [];
	scaninfo.trat       = [];
	scaninfo.imfm  	    = [dcminfo.Rows dcminfo.Columns];
	scaninfo.cntx       = [];
	scaninfo.cal        = [];
	scaninfo.min        = [];
	scaninfo.mag        = [];
	if isfield(dcminfo,'PixelSpacing')
        scaninfo.pixsize    = dcminfo.PixelSpacing';
        if length(scaninfo.pixsize) == 2; % pixel size info only for 2D 
           if isfield(dcminfo,'SliceThickness')
                scaninfo.pixsize(3) = dcminfo.SliceThickness;
           else
                scaninfo.pixsize(3) = 0;
           end
       elseif length(scaninfo.pixsize) < 2; % no pixel size info at all
           button = questdlg('No pixel size information was found! Do you want to continue to define them?',...
		'Continue Operation','Yes','No','No');
			if strcmp(button,'No')
                imaVOL = [];
                scaninfo = [];
                return;
            else
                prompt = {['Enter the X size in [mm] :'],['Enter the Y size in [mm] :'],['Enter the Z size in [mm] :']};
				dlg_title = 'Input for pixel sizes';
				num_lines= 1;
				def     = {'1','1','1'};
				scaninfo.pixsize  = str2double(inputdlg(prompt,dlg_title,num_lines,def))';
            end    
       elseif scaninfo.pixsize(3) == 0 % the Z size = 0
             button = questdlg('The Z pixel size is 0! Do you want to continue?',...
		'Continue Operation','Yes','No','No');
			if strcmp(button,'No')
                imaVOL = [];
                scaninfo = [];
                return;
            else
                prompt = {['Enter Z pixel size(mm) /X,Y size are '...
                            ,num2str(scaninfo.pixsize(1)),',',num2str(scaninfo.pixsize(2)),'/ :']};
				dlg_title = 'Input for pixel size';
				num_lines= 1;
				def     = {'1'};
				scaninfo.pixsize(3)  = str2double(inputdlg(prompt,dlg_title,num_lines,def));
            end
		end 
	else
        scaninfo.pixsize = [1 1];
	end
	if isfield(dcminfo,'NumberOfFrames')
        scaninfo.Frames             = dcminfo.NumberOfFrames;
        if isfield(dcminfo,'FrameTime')
            FrameTime = dcminfo.FrameTime/1000;%[sec];
            scaninfo.tissue_ts          = [FrameTime:FrameTime:scaninfo.Frames*FrameTime];
            scaninfo.start_times        = [];
            scaninfo.frame_lengths      = FrameTime*ones(1,scaninfo.Frames);
        else  %this is the case of SPECT projections: Frame means Projection.
            FrameTime = 1;
            scaninfo.tissue_ts          = [FrameTime:FrameTime:scaninfo.Frames*FrameTime];
            scaninfo.start_times        = [];
            scaninfo.frame_lengths      = FrameTime*ones(1,scaninfo.Frames);
        end
    elseif isfield(dcminfo,'NumberOfTimeSlices') %Philips Gemini TF dynamic dcm format
        scaninfo.Frames    = dcminfo.NumberOfTimeSlices;
        if scaninfo.Frames >1 
            scaninfo.tissue_ts          = zeros(1,scaninfo.Frames);
            scaninfo.start_times        = [];
            scaninfo.frame_lengths      = zeros(1,scaninfo.Frames);
            PhilipsDcmDyn = 1;
        end 
    else
        scaninfo.Frames =1;
        scaninfo.tissue_ts = [];
        scaninfo.start_times        = [];
        scaninfo.frame_lengths      = [];
	end
	if (num_of_files > 1) && (PhilipsDcmDyn == 0); 
        scaninfo.num_of_slice = num_of_files;
	elseif isfield(dcminfo,'NumberOfSlices')
        scaninfo.num_of_slice       = dcminfo.NumberOfSlices;
	else
        scaninfo.num_of_slice       = 1;
	end
	if scaninfo.Frames == scaninfo.num_of_slice
        scaninfo.Frames =1;
    end
    scaninfo.Frames = double(scaninfo.Frames);% for mia_gui FrameSlider
    scaninfo.num_of_slice = double(scaninfo.num_of_slice);% for mia_gui SliceSlider
	scaninfo.float = 0;
	scaninfo.FileType    = 'dcm';
	%
	% creating the imaVOL
	%
    
        
	if num_of_files > 1 
        delete(hm);
        if dcminfo.Rows <= 256 && isfield(dcminfo,'RescaleSlope') && ((dcminfo.RescaleSlope~=1) | (dcminfo.RescaleIntercept~=0)) ...
                && ( (mia_dicom_scaled == 1) | (mia_dicom_scaled == 2) )
            imaVOL = (zeros(dcminfo.Rows ,dcminfo.Columns,num_of_files));
        else
            imaVOL = (zeros(dcminfo.Rows ,dcminfo.Columns,num_of_files,'int16'));
        end
        % setup the progress bar
        info.color=[1 0 0];
		info.title='Creating image volume';
		info.size=1;
        info.pos='topleft';
		p=progbar(info);
		progbar(p,0);
        scaninfo.ImgIdx2FileIdx = zeros(num_of_files,1);
        for i=1:num_of_files
            filename = [dir_path,char(FileNames(i))];
            dcminfo = dicominfo(filename);
            if  ~isfield(dcminfo,'ManufacturerModelName')
                dcminfo.ManufacturerModelName = 'mia_noname';
            end
            if  ~isfield(dcminfo,'Units')
                dcminfo.Units = 'mia_noname';
            end
            SliceImage = dicomread(dcminfo);
            if size(SliceImage,3) > 1 % the slice image should be 2D
                % This is not true in the case of screencaptured dcm images
                disp('');
                disp('The dicom slices are 3D images!. They might be screencapture images.');  
                if ishandle(hm)
                    delete(hm);
                end
                imaVOL = []; scaninfo = []; dcminfo = [];
            end
            if isfield(dcminfo,'XrayTubeCurrent')
                scaninfo.mA(i) =  dcminfo.XrayTubeCurrent;
            end
            if isfield(dcminfo,'Exposure')
                scaninfo.exposure(i) =  dcminfo.Exposure;
            end
            if isfield(dcminfo,'CTDIvol')
                scaninfo.CTDIvol(i) =  dcminfo.CTDIvol;
            end
            if isfield(dcminfo,'EstimatedDoseSaving')
                scaninfo.EstimatedDoseSaving(i) =  dcminfo.EstimatedDoseSaving;
            end
            if isfield(dcminfo,'CTExposureSequence.Item_1.XrayTubeCurrentInmA')
                scaninfo.EstimatedDoseSaving(i) =  dcminfo.CTExposureSequence.Item_1.XrayTubeCurrentInmA;
            end
            if isfield(dcminfo,'RescaleSlope') && ...
                    (isfield(dcminfo,'ImageType') && isempty(strfind(dcminfo.ImageType,'LOCALIZER')))

                scaninfo.float = 1;
                if isfield(dcminfo,'Private_7053_1009') %Philips Gemini TF dcm SUV scaled format 
                        InjectActivity = dcminfo ...
                            .RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose; % [Bq]
                        InjectTimeString = dcminfo ...
                            .RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartTime;%[HHMMSS] 
                        if ~isempty(dcminfo.AcquisitionTime)
                            AcquisitionStartTimeStr = dcminfo.AcquisitionTime; %[HHMMSS] 
                        else
                            disp(' ');
                            disp(['Warning: dcminfo.AcquisitionTime is empty in the dcmInstance of ',num2str(dcminfo.InstanceNumber),' !']);
                            disp(['          The SUV calculation may not be correct!']);
                            AcquisitionStartTimeStr = dcminfo.SeriesTime;  %[HHMMSS]
                        end
                        AcquisitionStartDateStr = dcminfo.AcquisitionDate; %[YYYYMMSS]
                        ActualMidFrameTime = dcminfo.FrameReferenceTime/1000; %[sec]
                        ActivityTimeDuration = (datenum([AcquisitionStartDateStr,AcquisitionStartTimeStr],'yyyymmddHHMMSS') - ...
                            datenum([AcquisitionStartDateStr,InjectTimeString],'yyyymmddHHMMSS'))*24*3600 ...
                            + ActualMidFrameTime; %[sec] 
                        IsotopeHalfTime = dcminfo ...
                            .RadiopharmaceuticalInformationSequence.Item_1.RadionuclideHalfLife;%[sec]
                        if dcminfo.RescaleSlope == 1
                            BqCal = dcminfo.Private_7053_1009;
                        else
                            BqCal = 1; % ha a RescaleSlope~=1 (kb 2009 szep �ta), akkor a RescaleSlope is tartalmazza a 
                                        % Private_7053_1009 �rt�ket, teh�t az imaVOL sz�m�t�sakor m�r  elegend� 
                                        % a RescaleSlope haszn�lata.  
                        end
                        ActivityAtMidFrameTime =  InjectActivity*2^(-ActivityTimeDuration/IsotopeHalfTime);% [Bq]
                        %ActivityAtMidFrameTime =  InjectActivity;% [Bq] %for tests
                        Patientweight = dcminfo.PatientWeight*1000;%[g]
                        if mia_dicom_scaled == 1
                            SUVCal = BqCal/(ActivityAtMidFrameTime/Patientweight);
                        elseif mia_dicom_scaled == 2
                            SUVCal = BqCal; 
                        end
                elseif strcmp(dcminfo.ManufacturerModelName,'Gemini TF(C)') && ~isfield(dcminfo,'Private_7053_1009')
                        % Philips Gemini TF dcm non SUV scaled format. The
                        % pixeles contain the total counts, in that case
                        % the SUVcal factor make the counts -> cpm
                        % conversion
                        SUVCal = 1/(dcminfo.ActualFrameDuration/1000);
                elseif strcmp(dcminfo.ManufacturerModelName,'Discovery ST') && strcmp(dcminfo.Units,'BQML') ...
                        && isfield(dcminfo.RadiopharmaceuticalInformationSequence.Item_1,'RadionuclideTotalDose')
                        InjectActivity = dcminfo ...
                            .RadiopharmaceuticalInformationSequence.Item_1.RadionuclideTotalDose; % [Bq]
                        InjectTimeString = dcminfo ...
                            .RadiopharmaceuticalInformationSequence.Item_1.RadiopharmaceuticalStartTime;%[HHMMSS] 
                        AcquisitionStartTimeStr = dcminfo.SeriesTime; %[HHMMSS] 
                        AcquisitionStartDateStr = dcminfo.AcquisitionDate; %[YYYYMMSS]
                        ActivityTimeDuration = (datenum([AcquisitionStartDateStr,AcquisitionStartTimeStr],'yyyymmddHHMMSS') - ...
                            datenum([AcquisitionStartDateStr,InjectTimeString],'yyyymmddHHMMSS'))*24*3600; %[sec] 
                        % ActivityTimeDuration: dcminfo.DecayCorrection =
                        % 'START' a GE -n�l
                        IsotopeHalfTime = dcminfo ...
                            .RadiopharmaceuticalInformationSequence.Item_1.RadionuclideHalfLife;%[sec] 
                        ActivityAtScanStartTime =  InjectActivity*2^(-ActivityTimeDuration/IsotopeHalfTime);% [Bq]
                        Patientweight = dcminfo.PatientWeight*1000;%[g]
                        SUVCal = 1/(ActivityAtScanStartTime/Patientweight);
                else
                        SUVCal = 1;
                end
                tt(i) = dcminfo.InstanceNumber;
                fn{i} = filename;
                if mia_dicom_scaled == 1 | mia_dicom_scaled == 2
                    imaVOL(:,:,num_of_files + 3 - dcminfo.InstanceNumber) = squeeze( ...
                            double(SliceImage)*dcminfo.RescaleSlope + dcminfo.RescaleIntercept ...
                            )*SUVCal;
                else
                    
                    imaVOL(:,:,num_of_files + 3 - dcminfo.InstanceNumber) = squeeze(SliceImage*dcminfo.RescaleSlope + dcminfo.RescaleIntercept);
                end
                scaninfo.ImgIdx2FileIdx(num_of_files + 3 - dcminfo.InstanceNumber) = i;
            else
            
                imaVOL(:,:,num_of_files + 1 - dcminfo.InstanceNumber) = squeeze(SliceImage);
                scaninfo.ImgIdx2FileIdx(num_of_files + 1 - dcminfo.InstanceNumber) = i;
            end
        
            if (PhilipsDcmDyn == 1) && (mod(dcminfo.InstanceNumber,scaninfo.num_of_slice) == 0)
                scaninfo.tissue_ts(scaninfo.Frames + 1 - dcminfo.InstanceNumber/scaninfo.num_of_slice) = ...
                    dcminfo.FrameReferenceTime/1000;
                scaninfo.frame_lengths(scaninfo.Frames + 1 - dcminfo.InstanceNumber/scaninfo.num_of_slice)  = ...
                    dcminfo.ActualFrameDuration/1000;
            end
            if mod(i,3) == 0
                progbar(p,round(i/num_of_files*100));drawnow;
            end
            
        end
        % if some instancenumber were missing in the opened files the size(imaVOL,3)
        % might be increased above of num_of_files. The zeroed slices need to
        % be remove
        indexestodelete = [];
        if size(imaVOL,3) > num_of_files
            for i=1:size(imaVOL,3)
                imaVOLtmp = imaVOL(:,:,i);
                if isempty(find(imaVOLtmp))
                    indexestodelete = [indexestodelete i];
                end
            end
        end
        imaVOL(:,:,indexestodelete) = [];
        close(p);
    else
        SliceImage = dicomread(dcminfo);
        if isfield(dcminfo,'RescaleSlope')
            if size(SliceImage,3) > 1 % the slice image should be 2D
                % This is not true in the case of screencaptured dcm images
                disp('');
                disp('The dicom slices are 3D images!. They might be screencapture images.');  
                if ishandle(hm)
                    delete(hm);
                end
                imaVOL = []; scaninfo = []; dcminfo = [];
            end
            imaVOL = squeeze(int16( ...
                    double(SliceImage)*dcminfo.RescaleSlope + dcminfo.RescaleIntercept));
        else
            imaVOL = squeeze(SliceImage);
        end
        delete(hm);
    end

catch %in case of any error
    ErrorOnDicomOpening = lasterr
    if ishandle(hm)
        delete(hm);
    end
    if ishandle(p)
        delete(p);
    end
    imaVOL = []; scaninfo = []; dcminfo = [];
end