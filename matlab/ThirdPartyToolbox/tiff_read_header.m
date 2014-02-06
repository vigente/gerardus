function [TiffInfo,Img]=tiff_read_header(filename)
% This function tiff_read_header can be used to read the header(s) of
% complex tiff files. It implements a part of the Tiff version 5 and 6
% specifications.
%
%
% 	Info = tiff_read_header( Filename )
%		or
% 	[Info,I] = tiff_read_header( Filename )
%
% inputs,
%   Filename : filename of a Tif/Tiff file
%
% ouputs,
%   Info : Struct with Tiff tags, or in case of multiple images,
%          cell array with struct for every image.
%	I : Image Data (Only supported for non-compressed data)
%
% Note:
%   SubHeaders can be found in SubIFD1,SubIFD2 ...
%
%
% Literature,
%   http://www.awaresystems.be/imaging/tiff/tifftags/baseline.html
%   http://www.fileformat.info/format/tiff/corion.htm
%   http://www.compix.com/fileformattif.htm
%   http://partners.adobe.com/public/developer/en/tiff/TIFF6.pdf
%
% Written by D.Kroon 31-05-2012
%

% Create the Tiff-tag Dictionary
i=1;
Dic{i}='254,NewSubfileType'; i=i+1;
Dic{i}='255,SubfileType'; i=i+1;
Dic{i}='256,ImageWidth'; i=i+1;
Dic{i}='257,ImageHeight'; i=i+1;
Dic{i}='258,BitsPerSample'; i=i+1;
Dic{i}='259,Compression'; i=i+1;
Dic{i}='262,PhotometricInterpretation'; i=i+1;
Dic{i}='263,Treshholding'; i=i+1;
Dic{i}='264,CellWidth'; i=i+1;
Dic{i}='265,CellLength'; i=i+1;
Dic{i}='266,FillOrder'; i=i+1;
Dic{i}='269,DocumentName'; i=i+1;
Dic{i}='270,ImageDescription'; i=i+1;
Dic{i}='271,Make'; i=i+1;
Dic{i}='272,Model'; i=i+1;
Dic{i}='273,StripOffsets'; i=i+1;
Dic{i}='274,Orientation'; i=i+1;
Dic{i}='277,SamplesPerPixel'; i=i+1;
Dic{i}='278,RowsPerStrip'; i=i+1;
Dic{i}='279,StripByteCount'; i=i+1;
Dic{i}='280,MinSampleValue'; i=i+1;
Dic{i}='281,MaxSampleValue'; i=i+1;
Dic{i}='282,Xresolution'; i=i+1;
Dic{i}='283,Yresolution'; i=i+1;
Dic{i}='284,PlanarConfiguration'; i=i+1;
Dic{i}='285,PageName'; i=i+1;
Dic{i}='286,XPosition'; i=i+1;
Dic{i}='287,YPosition'; i=i+1;
Dic{i}='288,FreeOffsets'; i=i+1;
Dic{i}='289,FreeByteCounts'; i=i+1;
Dic{i}='290,GrayResponseUnit'; i=i+1;
Dic{i}='291,GrayResponseCurve'; i=i+1;
Dic{i}='292,T4Options'; i=i+1;
Dic{i}='293,T6Options'; i=i+1;
Dic{i}='296,ResolutionUnit'; i=i+1;
Dic{i}='297,PageNumber'; i=i+1;
Dic{i}='301,TransferFunction'; i=i+1;
Dic{i}='305,Software'; i=i+1;
Dic{i}='306,DateTime'; i=i+1;
Dic{i}='315,Artist'; i=i+1;
Dic{i}='316,HostComputer'; i=i+1;
Dic{i}='317,Predictor'; i=i+1;
Dic{i}='318,ColorImageType'; i=i+1;
Dic{i}='319,ColorList'; i=i+1;
Dic{i}='320,Colormap'; i=i+1;
Dic{i}='321,HalftoneHints'; i=i+1;
Dic{i}='322,TileWidth'; i=i+1;
Dic{i}='323,TileLength'; i=i+1;
Dic{i}='324,TileOffsets'; i=i+1;
Dic{i}='325,TileByteCounts'; i=i+1;
Dic{i}='326,BadFaxLines'; i=i+1;
Dic{i}='330,SubIFDs'; i=i+1;
Dic{i}='332,InkSet'; i=i+1;
Dic{i}='333,InkNames'; i=i+1;
Dic{i}='334,NumberOfInks'; i=i+1;
Dic{i}='336,DotRange'; i=i+1;
Dic{i}='337,TargetPrinter'; i=i+1;
Dic{i}='338,ExtraSamples'; i=i+1;
Dic{i}='339,SampleFormat'; i=i+1;
Dic{i}='340,SMinSampleValue'; i=i+1;
Dic{i}='341,MaxSampleValue'; i=i+1;
Dic{i}='342,TransferRange'; i=i+1;
Dic{i}='343,ClipPath'; i=i+1;
Dic{i}='33432,Copyright'; i=i+1;

% Split the Dictonary in Tag-ID and Tag-Name-list
TagId=zeros(1,length(Dic));
TagName=cell(1,length(Dic));
for j=1:length(Dic)
    bytes=uint8(Dic{j});
    bytes(bytes==32)=[];
    n=find(bytes==44);
    TagId(j)=str2double(char(bytes(1:n(1)-1)));
    TagName{j}=char(bytes(n(1)+1:end));
end

% Open the File
fid = fopen(filename,'r','l');

% Get File Size
fseek(fid,0,'eof');
fs = ftell(fid);


% Get the ByteOrder
fseek(fid,0,'bof');
ByteOrder= fread(fid,2,'char=>char')';
fclose(fid);
switch(ByteOrder)
    case 'II'
        fid = fopen(filename,'r','l');
    case 'MM'
        fid = fopen(filename,'r','b');
    otherwise
        fid = fopen(filename,'r','l');
end
fseek(fid,2,'bof');
VersionNumber= fread(fid,1,'uint16=>uint16');

% Get Position of the first Header (image file directory)
FirstIFD = fread(fid,1,'uint32=>uint32');

npage=1;
TiffInfo=cell(1,1);
TInfo=readIFD(fid,FirstIFD,TagId,TagName);
% If SubIFDs exist then there are sub images/headers
if(isfield(TInfo,'SubIFDs'))
    subIFD=TInfo.SubIFDs;
    for j=1:length(subIFD)
        TInfo.( ['SubIFD' num2str(j)])=readIFD(fid,subIFD(j),TagId,TagName);
    end
end
TiffInfo{npage}=TInfo;

% If NextIFD is larger then zero then there are multiple images/headers
NextIFD=TiffInfo{npage}.NextIFD;
while(NextIFD~=0)
    npage=npage+1;
    TInfo=readIFD(fid,NextIFD,TagId,TagName);
    % If SubIFDs exist then there are sub images/headers
    if(isfield(TInfo,'SubIFDs'))
        subIFD=TInfo.SubIFDs;
        for j=1:length(subIFD)
            TInfo.( ['SubIFD' num2str(j)])=readIFD(fid,subIFD(j),TagId,TagName);
        end
    end
    TiffInfo{npage}=TInfo;
    NextIFD=TiffInfo{npage}.NextIFD;
end
fclose(fid);


if(nargout>1)
    % Only works for No Compression
    if(TiffInfo{1}.Compression==1)
        fid = fopen(filename);
        % Read the Image Data
        for k=1:npage
            I=readImageData(fid,TiffInfo{k});
            if(k==1)
                Img=zeros([size(I) npage],class(I));
            end
            switch(ndims(I))
                case 2
                    Img(:,:,k)=I;
                case 3
                    Img(:,:,:,k)=I;
            end
            
        end
        fclose(fid)
    else
        Img=[];
    end
end

% If only one header, we output a struct instead of a cell
if(npage==1)
    TiffInfo=TiffInfo{npage};
end

% Read one Image File Directory (Header)
function TiffInfo=readIFD(fid,PositionIFD,TagId,TagName);
% Get number of Tags
fseek(fid,double(PositionIFD),'bof');
IfdLength=fread(fid,1,'uint16=>uint16');

TiffInfo=struct();
% Read all Tags
for j=1:IfdLength
    TagCode=fread(fid,1,'uint16=>uint16');
    TagType=fread(fid,1,'uint16=>uint16')';
    TagLength=fread(fid,1,'uint32=>uint32');
    
    switch(TagType)
        case 1 %byte
            nbyte=TagLength*1;
        case 2 % ASCII
            nbyte=TagLength*1;
        case 3 % Word
            nbyte=TagLength*2;
        case 4 % DWord - Uword
            nbyte=TagLength*4;
        case 5 % Rational (2 dwords, numerator and denominator)
            nbyte=TagLength*2*4;
        case 6 % 8-bit signed (twos-complement) integer.
            nbyte=TagLength;
        case 7 % A 8-bit byte undefined
            nbyte=TagLength*1;
        case 8 % 16-bit (2-byte) signed (twos-complement) integer.
            nbyte=TagLength*2;
        case 9 % A 32 bit(4-byte) signed (twos-complement) integer
            nbyte=TagLength*4;
        case 10 % Two SLONGs, numerator/denominator
            nbyte=TagLength*8;
        case 11 % Single precision (4-byte) IEEE format
            nbyte=TagLength*4;
        case 12 % Double precision (8-byte) IEEE format
            nbyte=TagLength*8;
        case 13 % uint32??
            nbyte=TagLength*4;
       otherwise

    end
    
    % If more bytes than 4, the data is stored
    % elsewhere in the file
    if(nbyte>4)
        TagDataOffset=fread(fid,1,'uint32=>uint32');
        cPos=ftell(fid);
        fseek(fid,double(TagDataOffset),'bof');
    end
    
    switch(TagType)
        case 1 %byte
            TagValue=fread(fid,TagLength,'uint8=>uint8');
        case 2 % ASCII
            TagValue=fread(fid,TagLength,'uint8=>uint8')';
            if(TagValue(end)==0), TagValue=TagValue(1:end-1); end
            TagValue=char(TagValue);
        case 3 % Word
            TagValue=fread(fid,TagLength,'uint16=>uint16');
        case 4 % DWord - Uword
            TagValue=fread(fid,TagLength,'uint32=>uint32');
        case 5 % Rational (2 dwords, numerator and denominator)
            TagValue=double(fread(fid,TagLength*2,'uint32=>uint32'));
            TagValue=TagValue(1:2:end)/TagValue(2:2:end);
        case 6 % An 8-bit (2-byte) signed (twos-complement) integer.
            TagValue=fread(fid,TagLength,'int8=>int8');
        case 7 % A 8-bit byte undefined
            TagValue=fread(fid,TagLength,'uint8=>uint8');
        case 8 % 16-bit (2-byte) signed (twos-complement) integer.
            TagValue=fread(fid,TagLength,'int16=>int16');
        case 9 % A 32 bit(4-byte) signed (twos-complement) integer
            TagValue=fread(fid,TagLength,'int32=>int32');
        case 10 % Two SLONGs, numerator/denominator
            TagValue=double(fread(fid,TagLength*2,'int32=>int32'));
            TagValue=TagValue(1:2:end)/TagValue(2:2:end);
        case 11 % Single precision (4-byte) IEEE format
            TagValue=fread(fid,TagLength,'single=>single');
        case 12 % Double precision (8-byte) IEEE format
            TagValue=fread(fid,TagLength,'double=>double');
        case 13
            TagValue=fread(fid,TagLength,'uint32=>uint32');

        otherwise
    end
    
    % If the data is less than 4 bytes it is zero padded
    % to 4 bytes
    if(nbyte<4)
        cPos=ftell(fid);
        fseek(fid,cPos+(4-nbyte),'bof');
    end
    
    % Go back from data position to tag positon
    if(nbyte>4)
        fseek(fid,cPos,'bof');
    end
    
    % Store Tag value in Struct will all tag-info
    n=find(TagId==TagCode);
    if(isempty(n))
        TiffInfo.( ['private_' num2str(TagCode)])=TagValue;
    else
        TName=TagName{n(1)};
        TiffInfo.(TName)=TagValue;
        
        % PhotometricInterpretation
        if(TagCode==262)
            switch(TagValue)
                case 0
                    TiffInfo.([TName 'String'])='GrayScaleWhite';
                case 1
                    TiffInfo.([TName 'String'])='GrayScaleBlack';
                case 2
                    TiffInfo.([TName 'String'])='RGB';
                case 3
                    TiffInfo.([TName 'String'])='PaletteColor';
                case 4
                    TiffInfo.([TName 'String'])='TransparencyMask';
                otherwise
                    TiffInfo.([TName 'String'])='Unknown';
            end
        end
        
        % Compression
        if(TagCode==259)
            switch(TagValue)
                case 1
                    TiffInfo.([TName 'String'])='NoCompression';
                case 2
                    TiffInfo.([TName 'String'])='Modified-Huffman-CCITT-Group3';
                case 3
                    TiffInfo.([TName 'String'])='Facsimile-compatible-CCITT-Group3';
                case 4
                    TiffInfo.([TName 'String'])='Facsimile-compatible-CCITT-Group4';
                case 5
                    TiffInfo.([TName 'String'])='LWZ';
                case 7
                    TiffInfo.([TName 'String'])='JPEG';
                case 8
                    TiffInfo.([TName 'String'])='ZIP';
                case 32773
                    TiffInfo.([TName 'String'])='PackBits';
                otherwise
                    TiffInfo.([TName 'String'])='Unknown';
            end
        end
    end
end
TiffInfo.NextIFD=fread(fid,1,'uint32=>uint32');

function Img=readImageData(fid,TiffInfo)
switch(TiffInfo.BitsPerSample(1))
    case 8
        nbyte=1;
        type='uint8';
    case 16
        nbyte=2;
        type='int16';
    case 32
        nbyte=4;
        type='single';
end

if(TiffInfo.PlanarConfiguration==1)
    s=double([ TiffInfo.SamplesPerPixel TiffInfo.ImageWidth TiffInfo.ImageHeight]);
else
    s=double([TiffInfo.ImageWidth TiffInfo.ImageHeight TiffInfo.SamplesPerPixel]);
end

Img=zeros([prod(s),1],type);
index=0;
for i=1:length(TiffInfo.StripOffsets);
    fseek(fid, TiffInfo.StripOffsets(i),'bof');
    data=fread(fid, TiffInfo.StripByteCount(i)/nbyte,type);
    Img((index+1):(index+length(data)))=data;
    index=index+length(data);
end
Img=reshape(Img,s);
if(TiffInfo.PlanarConfiguration==1)
    Img=permute(Img,[3 2 1]);
else
end



