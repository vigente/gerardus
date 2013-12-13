%% CBC Install for OPTI Toolbox
% Copyright (C) 2012 Jonathan Currie (I2C2)

% This file will help you compile Coin-Or Branch and Cut for use with 
% MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012

% To recompile you will need to get / do the following:

% 0) Compile CLP as per the instructions in opti_CLP_Install.m

% 1) Get CBC
% CBC is available from https://projects.coin-or.org/Cbc. Download 
% the source.

% 2) Compile CBC
% The easiest way to compile CBC is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

% %% Visual Studio Builder Commands
% clear
% path = 'full path to CBC here'; %e.g. 'C:\Solvers\CBC'
% n = 1;
% % CBC
% sdir = [path '\src'];
% hdrs = {[path '\..\CoinUtils\src'], [path '\..\Cgl\src'], [path '\..\Clp\src'], [path '\..\Osi\src'], [path '\..\BuildTools\headers']};
% name = 'libcbc';
% opts = [];
% opts.exPP = {'CBC_BUILD','COIN_FAST_CODE','CLP_FAST_CODE','USE_CBCCONFIG','COIN_NO_TEST_DUPLICATE','_CRT_SECURE_NO_WARNINGS'};
% opts.exclude = {'unitTest.cpp','unitTestClp.cpp','CoinSolve.cpp','CbcGeneric.cpp','CbcBranchBase.cpp'};
% VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% % CGL
% sdir = [path '\..\Cgl\src'];
% hdrs = {[path '\..\CoinUtils\src'], [path '\..\Clp\src\'], [path '\..\Osi\src'], [path '\..\BuildTools\headers']};
% name = 'libcgl';
% opts = [];
% opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};          
% VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% %Write the Solution File
% VS_WriteSol(VSPRJ)
% %%

% Once complete, you will have a directory called Cbc/libcbc. Open the
% Visual Studio 2012 project file, then complete the following steps:
%   a) Build a Win32 or x64 Release of each proejct to compile the code.
%   b) Copy the generated .lib files to the following folder:
%
%   OPTI/Solvers/cbc/Source/lib/win32 or win64
%
%   You will also need to copy all header files from Cbc/src to 
%   the following folder:
%
%   OPTI/Solvers/clp/Source/Include/Cbc
%
% 	And header files from Cgl/src (including sub-dirs) to:
%
%   OPTI/Solvers/clp/Source/Include/Cgl

% 3) CBC MEX Interface
% The CBC MEX Interface is a simple MEX interface I wrote to use CBC.

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the CBC MEX file. Once you have completed all 
% the above steps, simply run this file to compile CBC! You MUST BE in 
% the base directory of OPTI!

clear cbc

% Get Arch Dependent Library Path
libdir = opti_GetLibPath();
% Dependency Paths
clpdir = '..\..\clp\Source\';

fprintf('\n------------------------------------------------\n');
fprintf('CBC MEX FILE INSTALL\n\n');

%Get CBC Libraries
post = [' -IInclude\Cbc -IInclude\Osi -IInclude\Cgl -L' libdir ' -llibcbc -llibcgl -llibut'];
%Get CLP and Osi libraries
post = [post ' -I' clpdir 'Include -I' clpdir 'Include\Clp -I' clpdir 'Include\Coin -I' clpdir 'Include\Osi'];
post = [post ' -L' clpdir libdir ' -llibclp -llibcoinutils -llibosi -DCOIN_MSVS -output cbc'];

%CD to Source Directory
cdir = cd;
cd 'Solvers/cbc/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims cbcmex.cpp';
try
    eval([pre post])
    movefile(['cbc.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:clp','Error Compiling CBC!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');
