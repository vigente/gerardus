%% CoinUtils Install for OPTI Toolbox
% Copyright (C) 2012 Jonathan Currie (I2C2)

% This file will help you compile CoinUtils for use with MATLAB.

% My build platform:
% - Windows 7 SP1 x64
% - Visual Studio 2010

% The supplied MEX files will require the VC++ 2010 Runtime.

% You may download your own version of CoinUtils from 
% https://projects.coin-or.org/CoinUtils or I suggest you follow the
% instructions in opti_CLP_Install.m to download it as part of CLP. Note
% you will have to do the following changes to the CLP projects in order to
% be compatible with the code below, and the MEX interface:

% 1) Compile CoinUtils
% You must add the preprocessor "COIN_HAS_GLPK" to the C++ settings, and 
% add the include directory where the glpk.h include file is. This will
% require a download of GLPK. Then compile the libCoinUtils project.

% 2) Compile GLPK
% Follow the instructions in opti_GLPK_Install.m to build and compile
% GLPK. This is used for GMPL file reading.

% 3) MEX Interfaces
% The Read and Write MEX Interfaces are simple MEX interfaces I wrote to use 
% the CoinUtils File IO routines. They are supplied with the toolbox.

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the MEX files. Once you have completed all 
% the above steps, simply run this file to compile! You MUST BE in 
% the base directory of OPTI!

clear coinR coinW

% Get Arch Dependent Library Path
libdir = opti_GetLibPath();
% Dependency Paths
clpdir = '..\..\..\Solvers\clp\Source\';
glpkdir = '..\..\..\Solvers\glpk\Source\';

fprintf('\n------------------------------------------------\n');
fprintf('COINUTILS MEX FILE INSTALL\n\n');

%Get COIN Libraries (NOTE I have called CoinUtils with GMPL libCoinUtilsGMPL)!!
post = [' -I..\..\..\Solvers\clp\Source\Include\Coin -L' clpdir libdir ' -llibCoinUtilsGMPL'];
post = [post ' -L' glpkdir libdir ' -lglpk -DCOIN_MSVS -I' glpkdir '\Include\ '];

%CD to Source Directory
cdir = cd;
cd 'Utilities/File IO/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims coinR.cpp';
try
    eval([pre post])
     movefile(['coinR.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:coin','Error Compiling COINUTILS Read!\n%s',ME.message);
end
%Compile & Move
pre = 'mex -v -largeArrayDims coinW.cpp';
try
    eval([pre post])
     movefile(['coinW.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:coin','Error Compiling COINUTILS Write!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');
