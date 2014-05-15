%% CLP Install for OPTI Toolbox
% Copyright (C) 2013 Jonathan Currie (I2C2)

% This file will help you compile Coin-Or Linear Programming for use with 
% MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012

% To recompile you will need to get / do the following:

% 1) Get CLP
% CLP is available from http://www.coin-or.org/projects/Clp.xml. Download 
% the source.

% 2) Compile CLP & COIN Utils
% The easiest way to compile CLP is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

% %% Visual Studio Builder Commands
% clear
% path = 'full path to CLP here'; %e.g. 'C:\Solvers\CLP'
% glpkpath = [];%'full path to GLPK here'; %e.g. 'C:\Solvers\GLPK' (or leave blank)
% n = 1;
% %CLP
% sdir = [path '\src'];
% hdrs = {[path '\..\CoinUtils\src'], [path '\..\Osi\src'], [path '\..\BuildTools\headers']}; %[mumpspath '\libseq'], [mumpspath '\include'], 
% name = 'libclp';
% opts = [];
% opts.exPP = {'CLP_BUILD','_CRT_SECURE_NO_WARNINGS'};
% opts.exclude = {'MyEventHandler.cpp','MyMessageHandler.cpp','unitTest.cpp','CbcOrClpParam.cpp',...
%                 'ClpCholeskyMumps.cpp','ClpCholeskyUfl.cpp','ClpCholeskyWssmp.cpp','ClpCholeskyWssmpKKT.cpp','ClpMain.cpp'};
% opts.exFilter = {'Abc*','CoinAbc*'};
% VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% %CLP with Aboca (requires Intel C++ Compiler)
% sdir = [path '\src'];
% hdrs = {[path '\..\CoinUtils\src'], [path '\..\Osi\src'], [path '\..\BuildTools\headers']}; %[mumpspath '\libseq'], [mumpspath '\include'], 
% name = 'libclpabc';
% opts = [];
% opts.exPP = {'CLP_BUILD','_CRT_SECURE_NO_WARNINGS','CLP_HAS_ABC=4','__BYTE_ORDER=__LITTLE_ENDIAN','INTEL_COMPILER'};
% opts.exclude = {'MyEventHandler.cpp','MyMessageHandler.cpp','unitTest.cpp','CbcOrClpParam.cpp',...
%                 'ClpCholeskyMumps.cpp','ClpCholeskyUfl.cpp','ClpCholeskyWssmp.cpp','ClpCholeskyWssmpKKT.cpp','ClpMain.cpp'};
% VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% % CoinUtils
% sdir = [path '\..\CoinUtils\src'];
% hdrs = {[path '\..\BuildTools\headers']};
% name = 'libcoinutils';
% opts = [];
% opts.exPP = {'COINUTILS_BUILD','_CRT_SECURE_NO_WARNINGS'};          
% VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% % CoinUtils with GLPK (for GMPL CoinR reading)
% if(~isempty(glpkpath))
%     sdir = [path '\..\CoinUtils\src'];
%     hdrs = {[glpkpath '\src'], [path '\..\BuildTools\headers']};
%     name = 'libcoinutilsgmpl';
%     opts = [];
%     opts.exPP = {'COIN_HAS_GLPK','COINUTILS_BUILD','_CRT_SECURE_NO_WARNINGS'};          
%     VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% end
% % Osi
% sdir = [path '\..\Osi\src'];
% hdrs = {[path '\..\CoinUtils\src'] [path '\..\BuildTools\headers']};
% name = 'libosi';
% opts = [];
% opts.exPP = {'OSI_BUILD','_CRT_SECURE_NO_WARNINGS'};   
% opts.exFolder = {'OsiCpx','OsiGlpk','OsiGrb','OsiMsk','OsiSpx','OsiXpr','OsiCommonTest'};
% VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% %Write the Solution File
% VS_WriteSol(VSPRJ)
% %%

% Once complete, you will have a directory called Clp/libclp. Open the
% Visual Studio 2012 solution file, then complete the following steps:
%   a) John has made several bug fixes for the latest release, changesets
%   1963-1964 should all be manually included (or trunk version downloaded)
%   b) Build a Win32 or x64 Release of each project to compile the code.
%   c) Copy the generated .lib files to the following folder:
%
%   OPTI/Solvers/clp/Source/lib/win32 or win64
%
%   You will also need to copy all header files from Clp/src AND 
%   Clp/src/OsiClp to the following folder:
%
%   OPTI/Solvers/clp/Source/Include/Clp
%
% 	And header files from CoinUtils/src to:
%
%   OPTI/Solvers/clp/Source/Include/Coin
%
% 	And header files from Osi/src/Osi to:
%
%   OPTI/Solvers/clp/Source/Include/Osi

% 3) CLP MEX Interface
% The CLP MEX Interface is a simple MEX interface I wrote to use CLP.

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the CLP MEX file. Once you have completed all 
% the above steps, simply run this file to compile CLP! You MUST BE in 
% the base directory of OPTI!

clear clp

%Enable for Aboca Build (must compile using Intel C++ & VS2012 Linker)
haveABC = true;

% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('CLP MEX FILE INSTALL\n\n');

%Get CLP Libraries
post = [' -IInclude -IInclude/Clp -IInclude/Coin -L' libdir ' -llibCoinUtils -llibut -llibosi -DCOIN_MSVS -output clp'];
%Get Optional Aboca
if(haveABC)
    post = [post ' -llibclpabc -DINTEL_COMPILER -DCLP_HAS_ABC=4 -D__BYTE_ORDER=__LITTLE_ENDIAN']; 
else
    post = [post ' -llibclp'];
end

%CD to Source Directory
cdir = cd;
cd 'Solvers/clp/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims clpmex.cpp';
try
    eval([pre post])
     movefile(['clp.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:clp','Error Compiling CLP!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');
