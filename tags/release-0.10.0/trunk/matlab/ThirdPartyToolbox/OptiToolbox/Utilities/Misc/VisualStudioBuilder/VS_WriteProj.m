function [projPath,guid] = VS_WriteProj(srcpath,projName,incpath,opts)
%WRITEVS  Create a Visual Studio Project from selected paths
%
%   This function attempts to automatically create a Visual Studio Project
%   from a supplied source directory path (or paths). It works for me, but
%   may not for every project! Use with caution. Project will be created
%   one directory up from srcpath.
%
%   proj_path = VS_WriteProj(source_path,proj_name,include_path,options)
%
%   source_path     - absolute path to the source code to compile (may be a
%                     cell array of locations).
%   proj_name       - project name (optional).
%   include_path    - absolute path of directory to add to C/C++ Inlude
%                     Directories (may be a cell array of locations). Note
%                     files in these directories are not added to the
%                     project tree. (optional)
%   options         - structue with fields: (optional)
%                       'cpp'      read C++ and C files (true), read FORTRAN files (false) {true}
%                       'exPP'     cell array of extra preprocessor definitions {[]}
%                       'openMP'   true / false {false}
%                       'mkl'      true / false to include MKL headers {false}
%                       'charset'  Character Set {'unicode'}, 'multibyte'
%                       'toolset'  'v100' for VS2010, 'v110' for VS2012 {v110}
%                       'ifortver' Intel Fortran version {11.0}
%                       'exclude'  cell array of source files to exclude from project {[]}
%                       'exFilter' cell array of filtered source files (e.g. abc*) to exclude from project {[]}
%                       'exFolder' cell array of source folders to exclude from project {[]}

%Extract first source path if multiple paths specified (first is assumed base)
allpaths = srcpath;
if(iscell(srcpath))
    if(size(allpaths,2) > 1)
        allpaths = allpaths';
    end
    srcpath = allpaths{1};
end

%Error checking + default args
if(isempty(srcpath) || ~exist(srcpath,'dir')), error('Cannot find the directory:\n%s',srcpath); end
if(nargin < 2 || isempty(projName)), projName = '.lib'; end
if(nargin < 3), incpath = []; end
if(~isfield(opts,'cpp')), opts.cpp = true; end
if(~isfield(opts,'exPP')), opts.exPP = []; end
if(~isfield(opts,'openMP')), opts.openMP = false; end
if(~isfield(opts,'charset')), opts.charset = 'unicode'; end
if(~isfield(opts,'mkl')), opts.mkl = false; end
if(~isfield(opts,'toolset')), opts.toolset = 'v110'; end
if(~isfield(opts,'ifortver')), opts.ifortver = '11.0'; end
if(~isfield(opts,'exclude')), opts.exclude = []; end
if(~isfield(opts,'exFolder')), opts.exFolder = []; end
if(~isfield(opts,'exFilter')), opts.exFilter = []; end

%If VS2010 (vs100), remove from options (not required)
if(strcmpi(opts.toolset,'v100')), opts.toolset = []; end

switch(lower(opts.charset))
    case 'unicode'
        charset = 'Unicode';
    case 'multibyte'
        charset = 'MultiByte';
    otherwise
        error('Unknown character set selected, valid options are Unicode or Multibyte');
end

%Do the same for extra include paths (paths not added to the project tree)
incpaths = incpath;
if(iscell(incpath))
    if(size(incpaths,2) > 1)
        incpaths = incpaths';
    end
end

%Remove trailing \
if(iscell(allpaths))
    for i = 1:size(allpaths,1)
        if(allpaths{i}(end) == '\')
            allpaths{i}(end) = [];
        end
    end
else
    if(srcpath(end) == '\')
        srcpath(end) = [];
        allpaths(end) = [];
    end    
end
if(iscell(incpaths))
    for i = 1:size(incpaths,1)
        if(incpaths{i}(end) == '\')
            incpaths{i}(end) = [];
        end
    end
elseif(~isempty(incpaths))
    if(incpaths(end) == '\')
        incpaths(end) = [];
    end    
end

%Create Project Directory
if(any(strfind(srcpath,':'))) %assume full path
    %Go up one folder by removing last folder line
    ind = strfind(srcpath,'\');
    ipath = srcpath; ipath(ind(end):end) = [];
    %Make the directory
    mkdir(ipath,projName);
    projPath = [ipath '\' projName];
else    
    mkdir(cd,projName);    
    projPath = [cd '\' projName];
end

% THE REMAINDER OF THIS FILE IMPLEMENTS A STATIC LIBRARY WITH 32BIT AND 64BIT CONFIGURATIONS. 
% IT MAY NOT WORK WITH OLDER (OR FUTURE) VERSIONS OF VISUAL STUDIO

%Get a Project GUID
guid = getProjGUID();

%Header
if(opts.cpp)
    docNode = com.mathworks.xml.XMLUtils.createDocument('Project');
    p = docNode.getDocumentElement;
    p.setAttribute('DefaultTargets','Build');
    p.setAttribute('ToolsVersion','4.0');
    p.setAttribute('xmlns','http://schemas.microsoft.com/developer/msbuild/2003');
    %Project Configuration
    pc = createSection(docNode,'ItemGroup','ProjectConfigurations');
    pc.appendChild(writeProjConfig(docNode,'Debug','Win32'));
    pc.appendChild(writeProjConfig(docNode,'Debug','x64'));
    pc.appendChild(writeProjConfig(docNode,'Release','Win32'));
    pc.appendChild(writeProjConfig(docNode,'Release','x64'));
    p.appendChild(pc);    
    %Globals
    pc = createSection(docNode,'PropertyGroup','Globals');
    addElemText(docNode,pc,'ProjectGuid',['{' guid '}']); 
    addElemText(docNode,pc,'Keyword','Win32Proj');
    addElemText(docNode,pc,'RootNamespace',projName);
    p.appendChild(pc);
    %Import
    p.appendChild(createImport(docNode,'$(VCTargetsPath)\Microsoft.Cpp.Default.props',[],[]));    
    %Config Setup
    p.appendChild(createConfig(docNode,'Debug','Win32','StaticLibrary','true',charset,[],opts.toolset));
    p.appendChild(createConfig(docNode,'Debug','x64','StaticLibrary','true',charset,[],opts.toolset));
    p.appendChild(createConfig(docNode,'Release','Win32','StaticLibrary','false',charset,'true',opts.toolset));
    p.appendChild(createConfig(docNode,'Release','x64','StaticLibrary','false',charset,'true',opts.toolset));

    %Imports
    p.appendChild(createImport(docNode,'$(VCTargetsPath)\Microsoft.Cpp.props',[],[]));
    pc = docNode.createElement('ImportGroup');
    pc.setAttribute('Label','ExtensionSettings');
    p.appendChild(pc);
    p.appendChild(createPropSheet(docNode,'Debug','Win32'));
    p.appendChild(createPropSheet(docNode,'Debug','x64'));
    p.appendChild(createPropSheet(docNode,'Release','Win32'));
    p.appendChild(createPropSheet(docNode,'Release','x64'));

    %Macros
    pc = docNode.createElement('PropertyGroup');
    pc.setAttribute('Label','UserMacros');
    p.appendChild(pc);
else
    docNode = com.mathworks.xml.XMLUtils.createDocument('VisualStudioProject');
    p = docNode.getDocumentElement;
    p.setAttribute('ProjectType','typeStaticLibrary');
    p.setAttribute('ProjectCreator','Intel Fortran');
    p.setAttribute('Keyword','StaticLibrary');
    p.setAttribute('Version','11.0');
    p.setAttribute('ProjectIdGuid',['{' guid '}']);
    pc = createSection(docNode,'Platforms');
    pc.appendChild(createSection(docNode,'Platform','Name','Win32'));
    pc.appendChild(createSection(docNode,'Platform','Name','x64'));
    p.appendChild(pc);
end

%Read in Source + Header Files from Main Path(s)
[src,hdr] = VS_BuildFileList(allpaths,opts.exFolder,opts.cpp);
%Remove files in exclude list
if(~isempty(opts.exclude))
    for i = 1:size(src,1)
        if(~isempty(src{i,2}))
            ind = ismember(src{i,2},opts.exclude);
            if(any(ind))
                src{i,2}(ind) = [];
            end
        end
    end
end
%Remove files in exclude filter list
if(~isempty(opts.exFilter))
    for i = 1:size(src,1)
        for j = 1:length(opts.exFilter)
            if(~isempty(src{i,2}))
                ind = regexp(src{i,2},opts.exFilter{j});
                lind = true(size(ind));
                for k = 1:length(ind)
                    if(ind{k} == 1)
                        lind(k) = false;
                    end
                end
                src{i,2} = src{i,2}(lind);
            end
        end
    end
end

%Read in Extra Header Files to Add
projhdr = hdr;
[~,hdrinc] = VS_BuildFileList(incpaths,opts.exFolder,opts.cpp);
hdr = [hdr;hdrinc];

% %Add MKL include path if requested
if(opts.mkl)
    [~,~,mkl_inc] = opti_FindMKL();
    if(isempty(hdr))
        hdr = {mkl_inc {'mkl.h'} 1};
    else
        no = hdr{end,end}+1;
        hdr = [hdr;{mkl_inc {'mkl.h'} no}];
    end
end

if(opts.cpp)
    %Debug Detailed Settings
    p.appendChild(writeDebugDetail(docNode,allpaths,hdr,'Win32',opts.exPP,opts.openMP));
    p.appendChild(writeDebugDetail(docNode,allpaths,hdr,'x64',opts.exPP,opts.openMP));
    %Release Detailed Settings
    p.appendChild(writeReleaseDetail(docNode,allpaths,hdr,'Win32',opts.exPP,opts.openMP));
    p.appendChild(writeReleaseDetail(docNode,allpaths,hdr,'x64',opts.exPP,opts.openMP));

    %Write VS Filters
    VS_WriteFilters(projPath,projName,allpaths,src,projhdr);

    %Source Files
    p.appendChild(createFileList(docNode,allpaths,src,'ClCompile'));
    %Header Files
    p.appendChild(createFileList(docNode,allpaths,hdr,'ClInclude'));

    %Imports
    p.appendChild(createImport(docNode,'$(VCTargetsPath)\Microsoft.Cpp.targets',[],[]));
    pc = docNode.createElement('ImportGroup');
    pc.setAttribute('Label','ExtensionTargets');
    p.appendChild(pc);
else
    %Configurations
    pc = createSection(docNode,'Configurations');
    pc.appendChild(writeIFortConfig(docNode,'Debug','Win32',opts.exPP,incpaths));
    pc.appendChild(writeIFortConfig(docNode,'Release','Win32',opts.exPP,incpaths));
    pc.appendChild(writeIFortConfig(docNode,'Debug','x64',opts.exPP,incpaths));
    pc.appendChild(writeIFortConfig(docNode,'Release','x64',opts.exPP,incpaths));
    p.appendChild(pc);
    %Write VS Filters
%     VS_WriteFilters(projPath,projName,allpaths,src,projhdr);
    %Files
    pc = createSection(docNode,'Files');
    pc.appendChild(createMulSection(docNode,'Filter',{{'Name','Header Files'},{'Filter','fi;fd'}}));
    pc.appendChild(createMulSection(docNode,'Filter',{{'Name','Resource Files'},{'Filter','rc;ico;cur;bmp;dlg;rc2;rct;bin;rgs;gif;jpg;jpeg;jpe'}}));
    pc.appendChild(createMulSection(docNode,'Filter',{{'Name','Source Files'},{'Filter','f90;for;f;fpp;ftn;def;odl;idl;F'}}));
    createFileList(docNode,allpaths,src,'File',true,pc);
%     createFileList(docNode,allpaths,hdr,'File',true,pc);
    p.appendChild(pc);
end

%Write Whole Project
xmlwrite([projPath '\' projName '.xml'],docNode);
if(opts.cpp)
    xmlwrite([projPath '\' projName '.vcxproj'],docNode);
else
    xmlwrite([projPath '\' projName '.vfproj'],docNode);
end

%Create File List Section (lists files to include / compile)
function pc = createFileList(docNode,upath,files,elem,isFort,pRoot)
if(nargin < 5), isFort = false; end
if(nargin < 6), pRoot = false; end
if(~isFort)
    pc = docNode.createElement('ItemGroup');
end
for i = 1:size(files,1) %for each path
    if(iscell(upath))
        p = upath{files{i,3}};
    else
        p = upath;
    end

    %Check we are on the same path as our base one
    h = files{i,1};
    if(isempty(strfind(h,p))) %empty means not same base bath
        attp = getPathSep(h,p);
    else %OK contained with base path
        if(any(strfind(p,':')))
            ind = strfind(p,'\'); %find end file separator
            len2 = ind(end)+1; %drop last folder (we are always back one)
        else
            len2 = 1;
        end    
        attp = files{i,1}(len2:end);
    end
    if(~isempty(files{i,2}))
        for j = 1:size(files{i,2},1) %for each file
            if(~isFort)
                pc1 = docNode.createElement(elem);
                pc1.setAttribute('Include',['..\' attp '\' files{i,2}{j}]);
                pc.appendChild(pc1);
            else
                pc = docNode.createElement(elem);
                pc.setAttribute('RelativePath',['..\' attp '\' files{i,2}{j}]);
                pRoot.appendChild(pc);
            end
        end
    end
end

%Create an XML Section
function pc = createSection(docNode,name,labelname,labelval)
pc = docNode.createElement(name);
switch(nargin)
    case 3
        pc.setAttribute('Label',labelname);
    case 4
        pc.setAttribute(labelname,labelval);
end

%Create an XML Section with multiple attributes
function pc = createMulSection(docNode,name,attr)
pc = docNode.createElement(name);
for i = 1:length(attr)
    pc.setAttribute(attr{i}{1},attr{i}{2});
end

%Lowlevel routine to add an element with text
function addElemText(docNode,pc,elem,text)
c = docNode.createElement(elem);
c.appendChild(docNode.createTextNode(text));
pc.appendChild(c);

%Create Build Configuration Section
function pc = createConfig(docNode,config,plat,lib,debug,char,opt,ts)
pc = docNode.createElement('PropertyGroup');
pc.setAttribute('Condition',['''$(Configuration)|$(Platform)''==''' config '|' plat '''']);
pc.setAttribute('Label','Configuration');
if(~isempty(lib));   addElemText(docNode,pc,'ConfigurationType',lib); end
if(~isempty(debug)); addElemText(docNode,pc,'UseDebugLibraries',debug); end
if(~isempty(char));  addElemText(docNode,pc,'CharacterSet',char); end
if(~isempty(opt));   addElemText(docNode,pc,'WholeProgramOptimization',opt); end
if(~isempty(ts)),    addElemText(docNode,pc,'PlatformToolset',ts); end

%Create Import Section
function pc = createImport(docNode,name,cond,label)
pc = docNode.createElement('Import');
pc.setAttribute('Project',name);
if(~isempty(cond)); pc.setAttribute('Condition',cond); end
if(~isempty(label)); pc.setAttribute('Label',label); end

%Create Property Sheet Section
function pc = createPropSheet(docNode,config,plat)
pc = docNode.createElement('ImportGroup');
pc.setAttribute('Label','PropertySheets');
pc.setAttribute('Condition',['''$(Configuration)|$(Platform)''==''' config '|' plat '''']);
pc.appendChild(createImport(docNode,'$(UserRootDir)\Microsoft.Cpp.$(Platform).user.props','exists(''$(UserRootDir)\Microsoft.Cpp.$(Platform).user.props'')','LocalAppDataPlatform'));

%Write Project Configuration Section
function pc = writeProjConfig(docNode,config,plat)
pc = docNode.createElement('ProjectConfiguration');
pc.setAttribute('Include',[config '|' plat]);
addElemText(docNode,pc,'Configuration',config)
addElemText(docNode,pc,'Platform',plat);

%Write Fortran Configuration Section
function pc = writeIFortConfig(docNode,config,plat,exPP,exInc)
pc = docNode.createElement('Configuration');
pc.setAttribute('Name',[config '|' plat]);
pc.setAttribute('ConfigurationType','typeStaticLibrary');
tl = docNode.createElement('Tool');
tl.setAttribute('Name','VFFortranCompilerTool');
tl.setAttribute('SuppressStartupBanner','true');
tl.setAttribute('Preprocess','preprocessYes');
tl.setAttribute('EnableEnhancedInstructionSet','codeArchSSE2');
if(nargin > 3 && ~isempty(exPP))
    tl.setAttribute('PreprocessorDefinitions',concatenatePP(exPP));
end	
if(nargin > 4 && ~isempty(exInc))
    tl.setAttribute('AdditionalIncludeDirectories',concatenatePP(exInc));
end
if(strcmpi(config,'debug'))
    tl.setAttribute('DebugInformationFormat','debugEnabled');
    tl.setAttribute('Optimization','optimizeDisabled');
    tl.setAttribute('RuntimeLibrary','rtMultiThreadedDebugDLL');
    tl.setAttribute('Traceback','true');
    tl.setAttribute('BoundsCheck','true');
    tl.setAttribute('StackFrameCheck','true');
else
    tl.setAttribute('Optimization','optimizeFull');
    tl.setAttribute('RuntimeLibrary','rtMultiThreadedDLL');
end
pc.appendChild(tl);
if(strcmpi(plat,'x64'))
    tl = docNode.createElement('Tool');
    tl.setAttribute('Name','VFMidlTool');
    tl.setAttribute('SuppressStartupBanner','true');
    tl.setAttribute('TargetEnvironment','midlTargetAMD64');
    pc.appendChild(tl);
end


%Write C++ Debug Detail Section
function pc = writeDebugDetail(docNode,upath,hdr,plat,exPP,openMP)
pc = docNode.createElement('ItemDefinitionGroup');
pc.setAttribute('Condition',['''$(Configuration)|$(Platform)''==''Debug|',plat,'''']);
cl = docNode.createElement('ClCompile');
cl.appendChild(docNode.createElement('PrecompiledHeader'));
addElemText(docNode,cl,'WarningLevel','Level3');
addElemText(docNode,cl,'Optimization','Disabled');
if(openMP), addElemText(docNode,cl,'OpenMPSupport','true'); end
if(~isempty(exPP)), exPP = concatenatePP(exPP); end
addElemText(docNode,cl,'PreprocessorDefinitions',['WIN32;_DEBUG;_LIB;' exPP '%(PreprocessorDefinitions)']);
addElemText(docNode,cl,'AdditionalIncludeDirectories',[includeStr(upath,hdr) '%(AdditionalIncludeDirectories)']);
pc.appendChild(cl);
lk = docNode.createElement('Link');
addElemText(docNode,lk,'SubSystem','Windows');
addElemText(docNode,lk,'GenerateDebugInformation','true');
pc.appendChild(lk);

%Write C++ Release Detail Section
function pc = writeReleaseDetail(docNode,upath,hdr,plat,exPP,openMP)
pc = docNode.createElement('ItemDefinitionGroup');
pc.setAttribute('Condition',['''$(Configuration)|$(Platform)''==''Release|',plat,'''']);
cl = docNode.createElement('ClCompile');
cl.appendChild(docNode.createElement('PrecompiledHeader'));
addElemText(docNode,cl,'WarningLevel','Level3');
addElemText(docNode,cl,'Optimization','MaxSpeed');
addElemText(docNode,cl,'IntrinsicFunctions','true');
if(~strcmpi(plat,'x64')), addElemText(docNode,cl,'EnableEnhancedInstructionSet','StreamingSIMDExtensions2'); end
addElemText(docNode,cl,'MultiProcessorCompilation','true');
addElemText(docNode,cl,'RuntimeTypeInfo','true');
addElemText(docNode,cl,'FavorSizeOrSpeed','Speed');
addElemText(docNode,cl,'InlineFunctionExpansion','OnlyExplicitInline');
addElemText(docNode,cl,'RuntimeLibrary','MultiThreadedDLL');
if(openMP), addElemText(docNode,cl,'OpenMPSupport','true'); end
if(~isempty(exPP)), exPP = concatenatePP(exPP); end
addElemText(docNode,cl,'PreprocessorDefinitions',['WIN32;NDEBUG;_LIB;' exPP '%(PreprocessorDefinitions)']);
addElemText(docNode,cl,'AdditionalIncludeDirectories',[includeStr(upath,hdr) '%(AdditionalIncludeDirectories)']);
pc.appendChild(cl);
lk = docNode.createElement('Link');
addElemText(docNode,lk,'SubSystem','Windows');
addElemText(docNode,lk,'GenerateDebugInformation','false');
addElemText(docNode,lk,'OptimizeReferences','true');
pc.appendChild(lk);

%Function to get a GUID from .EXE
function str = getProjGUID()
cdir = cd;
cd(GETCD);
dos(['"' GETCD 'genGUID" 1']);
cd(cdir);
fh = fopen([GETCD 'guids.txt']);
str = fgetl(fh);
if(length(str) < 5)
    error('Error reading GUID');
end
fclose(fh);

%Generate an include string from supplied cell array
function str = includeStr(upath,hdr)
str = [];
for i = 1:size(hdr,1) %for each path
    if(iscell(upath))
        p = upath{hdr{i,3}};
    else
        p = upath;
    end
    %Check we are on the same path as our base one
    h = hdr{i,1};
    if(isempty(strfind(h,p))) %empty means not same base path
        attp = getPathSep(h,p);
    else %OK contained with base path
        if(any(strfind(p,':')))
            ind = strfind(p,filesep); %find end file separator
            len2 = ind(end)+1; %drop last folder (we are always back one)
        else
            len2 = 1;
        end    
        attp = hdr{i,1}(len2:end);
    end
    
    str = [str '..' filesep attp ';']; %#ok<AGROW>
end

%Determine how many ../ we need and get path
function str = getPathSep(h,p)
%Run forward until we find path section that doesn't match
for i = 1:min(length(h),length(p))
    if(~strcmp(h(1:i),p(1:i)))
        break;
    end
end
%Go back until last pathsep
for k = i:-1:1
    if(strcmp(h(k),filesep))
        break;
    end
end
i = k+1;
%Now determine how many file seps in front of us
ind = strfind(p(i:end),filesep);
%Determine how many explicit steps back
indb = strfind(p(i:end),['..' filesep]);
str = '';
for j = 1:(length(ind)-length(indb)*2) %check the *2 (or mod of) with libClp
    str = sprintf('%s..%s',str,filesep);
end
%Concatenate with new bit 
str = [str h(i:end)];

%Get cwd of this file
function str = GETCD()
str = which('VS_WriteProj.m');
ind = strfind(str,filesep);
str(ind(end)+1:end) = [];

%Concatentate preprocessors
function str = concatenatePP(pp)
if(isempty(pp))
    str = [];
    return;
end
if(iscell(pp))
    str = [];
    for i = 1:length(pp)
        str = sprintf('%s%s;',str,pp{i});
    end
elseif(ischar(pp))
    str = pp;
else
    error('Unknown preprocessor format');
end
