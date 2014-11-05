function contourStruct = getCMR42Contours(xmlFile,dicomdirpath)
%GETCONTOURS   extracts pre-defined contours from a cmr42-generated xml
%   file, for a specified DICOM series. This function depends on Laszlo
%   Balkay's excellent DICOMDIR Reader
%   (http://www.mathworks.com/matlabcentral/fileexchange/7926-dicomdir-reader)
%   to select the series; the presence of an appropriate DICOMDIR file is
%   necessary for this to work.
%
%   Usage:
%     GETCONTOURS()
%     GETCONTOURS(xmlFile)
%     GETCONTOURS(xmlFile,dicomdirpath)
%
%   - xmlFile is a string corresponding to the cmr42-generated xml file
%   - dicomdirpath is a string corresponding to the path where a DICOMDIR is
%     present
%
%   If the function is called with no arguments, the user can select the
%   xmlFile and dicomdirpath using a file dialogue. Similarly, the function
%   can be called with the xmlFile argument only
%
%   The function returns a struct whose fields are any contour names found
%   in the xml file for the particular study in question
%
%   Example:
%     >> contourStruct = GETCONTOURS(xmlFile,dicomdirpath)
% 
%     contourStruct = 
% 
%              laxLvExtentPoints: [256x176x25 double]
%           saendocardialContour: [256x176x25 double]
%            saepicardialContour: [256x176x25 double]
%             sapapilMuscContour: [256x176x25 double]
%         sapapilMuscContour0001: [256x176x25 double]
%         sapapilMuscContour0002: [256x176x25 double]
%         sapapilMuscContour0003: [256x176x25 double]
%
%     >> imaVOL = loaddcm(loaddcmdir(dicomdirpath));
%     >> subimage(mat2gray(imaVOL(:,:,1))); hold on
%     >> imagesc(contourStruct.saepicardialContour(:,:,1),'AlphaData',contourStruct.saepicardialContour(:,:,1));
%
%   See also LOADDCM, LOADDCMDIR.

% Author: Tasos Papastylianou <tasos.papastylianou@gmail.com>
% Copyright © 2013-2014 University of Oxford
% Version: 0.1
% $Rev$
% $Date$
%
% University of Oxford means the Chancellor, Masters and Scholars of
% the University of Oxford, having an administrative office at
% Wellington Square, Oxford OX1 2JD, UK.
%
% This file is part of Gerardus.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. The offer of this
% program under the terms of the License is subject to the License
% being interpreted in accordance with English Law and subject to any
% action against the University of Oxford being under the jurisdiction
% of the English Courts.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% ======================== Note on DICOM terminology ======================
% Unfortunately, while the Study UID is present, the Series UID does not
% appear anywhere in the cmr42-generated xml file for any of the series
%
% In the context of DICOM, a study is a collection of series, it is the top
% of the hierarchy; while a DICOMDIR may contain information about many
% studies, usually a DICOMDIR describing a collection of DICOM images will
% deal with a single study, such as a series of cardiovascular
% investigations for a patient.)
%
% In the context of DICOM, a series is a collection of DICOM images,
% representing a single investigation. A typical example would be a
% collection of, say, 25 DICOM images, corresponding to 25 time frames of
% an SA cineMRI image at a particular level. If the entire heart is
% captured in 10 SA levels, there would usually be 10 'series' consisting
% of 25 images (i.e. frames) each.
%
% Each single DICOM image is identified by its filename, and also uniquely
% by its SOP Identifier; the cmr42-generated xml file only uses the latter
% for identification
% =========================================================================

  if nargin == 0
    [xmlFileName xmlPathName] = uigetfile('*.*','Select cmr42-generated XML file to Open');
    xmlFile = fullfile(xmlPathName,xmlFileName);
    [dcmSeries, dcmPatient, SeriesList, SeriesListNumInfo, ImageInfos, SeriesInfos] = loaddcmdir;
  elseif nargin == 1
    [dcmSeries, dcmPatient, SeriesList, SeriesListNumInfo, ImageInfos, SeriesInfos] = loaddcmdir;
  elseif nargin == 2
    [dcmSeries, dcmPatient, SeriesList, SeriesListNumInfo, ImageInfos, SeriesInfos] = loaddcmdir(dicomdirpath);
  end
    
  global dcmdirlistVal; % will have been generated in the base workspace by loaddcmdir
  dicomfiles = dcmPatient.Study.Series(dcmdirlistVal).ImageNames;
  dicomUIDs = dcmPatient.Study.Series(dcmdirlistVal).ReferencedSOPInstanceUIDInFile;
  studyUID = getStudyUID(dcmSeries);
  seriesLength = length(dicomfiles);

  xmlDoc = xmlread(xmlFile);
  xmlRoot = xmlDoc.getDocumentElement();
  xmlStudyMapStatesNode = getStudyMapStatesNode(xmlRoot);
  xmlStudyMapStatesNodeListChildren = getStudyMapStatesNodeListChildren(xmlStudyMapStatesNode);
  xmlCorrectStudyUidNode = getCorrectStudyUidNode(xmlStudyMapStatesNodeListChildren,studyUID);
  xmlCorrectImageStatesNode = getCorrectImageStatesNode(xmlCorrectStudyUidNode);
  
  contourNames = getContourNames(xmlCorrectImageStatesNode);
  contourVOL = zeros(getRows(dcmSeries), getColumns(dcmSeries),seriesLength);
  contourStruct = initialiseContourStruct(contourNames,contourVOL);
    
  for seriesFrame = 1 : seriesLength 
    xmlCorrectFrameNode = getCorrectFrameNode(xmlCorrectImageStatesNode,dicomUIDs{seriesFrame});
    for contourNameIndex = 1 : length(contourNames)
      contourName = contourNames{contourNameIndex};
      contourStruct.(contourName)(:,:,seriesFrame) = collectContourPoints(xmlCorrectFrameNode,contourVOL(:,:,1),contourName);
    end
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%
% HELPER FUNCTIONS BELOW
%%%%%%%%%%%%%%%%%%%%%%%%

function rows = getRows(dcmSeries)
  dcmSeriesFirstImage = dcmSeries;
  dcmSeriesFirstImage.Images = dcmSeries.Images(1);
  [imaVOL, scaninfo, dcminfo] = loaddcm(dcmSeriesFirstImage);
  rows = dcminfo.Rows;
end

function columns = getColumns(dcmSeries)
  dcmSeriesFirstImage = dcmSeries;
  dcmSeriesFirstImage.Images = dcmSeries.Images(1);
  [imaVOL, scaninfo, dcminfo] = loaddcm(dcmSeriesFirstImage);
  columns = dcminfo.Columns;
end

function studyUID = getStudyUID(dcmSeries)
  dcmSeriesFirstImage = dcmSeries;
  dcmSeriesFirstImage.Images = dcmSeries.Images(1);
  [imaVOL, scaninfo, dcminfo] = loaddcm(dcmSeriesFirstImage);
  studyUID = dcminfo.StudyInstanceUID;
end

function studyMapStatesNode = getStudyMapStatesNode(Root)
  for i = 0 : (Root.getChildNodes.getLength - 1)
    if Root.getChildNodes.item(i).hasAttributes
      if strcmp(Root.getChildNodes.item(i).getAttribute('Hash:key'),'StudyMapStates')
        studyMapStatesNode = Root.getChildNodes.item(i);
        return
      end
    end
  end
end

function studyMapStatesNodeListChildren = getStudyMapStatesNodeListChildren(studyMapStatesNode)
  studyMapStatesNodeListChildren = cell(0);
  for i = 0 : (studyMapStatesNode.getChildNodes.getLength - 1)
    if strcmp(studyMapStatesNode.getChildNodes.item(i).getNodeName,'List:item')
      studyMapStatesNodeListChildren{end+1} = studyMapStatesNode.getChildNodes.item(i);
    end
  end
end

function correctStudyUidNode = getCorrectStudyUidNode(studyMapStatesNodeListChildren,studyUID)
  for i = 1 : length(studyMapStatesNodeListChildren)
    for j = 0 : (studyMapStatesNodeListChildren{i}.getChildNodes.getLength - 1)
      if studyMapStatesNodeListChildren{i}.getChildNodes.item(j).hasAttributes
        if studyMapStatesNodeListChildren{i}.getChildNodes.item(j).hasAttribute('Hash:key')
          if strcmp(studyMapStatesNodeListChildren{i}.getChildNodes.item(j).getAttribute('Hash:key'),'StudyUid')
            if studyMapStatesNodeListChildren{i}.getChildNodes.item(j).getChildNodes.getLength == 1       
              if strcmp(studyMapStatesNodeListChildren{i}.getChildNodes.item(j).getChildNodes.item(0).getNodeName,'#text')
                if strcmp(studyMapStatesNodeListChildren{i}.getChildNodes.item(j).getChildNodes.item(0).getData,studyUID)
                  correctStudyUidNode = studyMapStatesNodeListChildren{i}.getChildNodes.item(j);
                  return
                end
              end
            end
          end
        end
      end
    end
  end
end

function correctImageStatesNode = getCorrectImageStatesNode(correctStudyUidNode)
  tempnode = correctStudyUidNode;
  while ~strcmp(class(tempnode),'double')
    tempnode = tempnode.getNextSibling;
    if strcmp(tempnode.getNodeName,'Hash:item')
      if tempnode.hasAttributes
        if strcmp(tempnode.getAttribute('Hash:key'),'ImageStates')
          correctImageStatesNode = tempnode;
          return
        end
      end
    end
  end
end

function correctFrameNode = getCorrectFrameNode(correctImageStatesNode,dicomUIDspecifier)
  for i = 0 : (correctImageStatesNode.getChildNodes.getLength - 1)
    if strcmp(correctImageStatesNode.getChildNodes.item(i).getNodeName,'Hash:item')
      if strcmp(correctImageStatesNode.getChildNodes.item(i).getAttribute('Hash:key'),dicomUIDspecifier)
        correctFrameNode = correctImageStatesNode.getChildNodes.item(i);
        return
      end
    end
  end
end

function contourNames = getContourNames(correctImageStatesNode)
  contourNames = cell(0);
  for i = 0 : (correctImageStatesNode.getChildNodes.getLength - 1)
    frameNode = correctImageStatesNode.getChildNodes.item(i);
    for j = 0 : (frameNode.getChildNodes.getLength - 1)
      if frameNode.getChildNodes.item(j).hasAttributes
        if strcmp(frameNode.getChildNodes.item(j).getAttribute('Hash:key'),'Contours')
          contourNode = frameNode.getChildNodes.item(j);
          if str2num(contourNode.getAttribute('Hash:count')) == 0
            continue
          else
            for k = 0 : (contourNode.getChildNodes.getLength - 1)
              if strcmp(contourNode.getChildNodes.item(k).getNodeName,'Hash:item')
                contourNames{end+1} = contourNode.getChildNodes.item(k).getAttribute('Hash:key');
                contourNames{end} = contourNames{end}.toCharArray;
                contourNames{end} = contourNames{end}';
              end
            end
          end
        end
      end
    end  
  end
  contourNames = unique(contourNames);
end

function contourStruct = initialiseContourStruct(contourNames,contourVOL)
  contourStruct = struct();
  for i = 1 : length(contourNames)
    contourStruct.(contourNames{i}) = contourVOL;
  end
end

function contourImg = collectContourPoints(correctFrameNode,contourImg,contourName)
  for i = 0 : (correctFrameNode.getChildNodes.getLength - 1)
    if correctFrameNode.getChildNodes.item(i).hasAttributes
      if strcmp(correctFrameNode.getChildNodes.item(i).getAttribute('Hash:key'),'Contours')
        contourNode = correctFrameNode.getChildNodes.item(i);
        if str2num(contourNode.getAttribute('Hash:count')) == 0
          return % returns unchanged the zero image that was supplied to the function
        else
          for j = 0 : ( contourNode.getChildNodes.getLength - 1 )
            if strcmp(contourNode.getChildNodes.item(j).getNodeName,'Hash:item')
              if strcmp(contourNode.getChildNodes.item(j).getAttribute('Hash:key'),contourName)
                contourImg = getPoints(contourImg,contourNode.getChildNodes.item(j));
              end
            end
          end
        end
      end
    end
  end
end

function contourImg = getPoints(contourImg,contourNameNode)
  for i = 0 : (contourNameNode.getChildNodes.getLength - 1)
    if strcmp(contourNameNode.getChildNodes.item(i).getNodeName,'Hash:item')
      if strcmp(contourNameNode.getChildNodes.item(i).getAttribute('Hash:key'),'Points')
        pointsNode = contourNameNode.getChildNodes.item(i);
      end
      if strcmp(contourNameNode.getChildNodes.item(i).getAttribute('Hash:key'),'SubpixelResolution')
        subPixelResolution = str2num(contourNameNode.getChildNodes.item(i).getTextContent);
      end      
    end
  end
  
  for i = 0 : ( pointsNode.getChildNodes.getLength - 1)
    if strcmp(pointsNode.getChildNodes.item(i).getNodeName,'List:item')
      listitem = pointsNode.getChildNodes.item(i);
      x = str2num(listitem.getElementsByTagName('Point:x').item(0).getTextContent)/subPixelResolution;
      y = str2num(listitem.getElementsByTagName('Point:y').item(0).getTextContent)/subPixelResolution;
      contourImg(round(y),round(x)) = 1;
    end
  end
end













