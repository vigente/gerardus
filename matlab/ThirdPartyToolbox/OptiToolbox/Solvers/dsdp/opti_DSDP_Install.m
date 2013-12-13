%% DSDP Install for OPTI Toolbox
% Copyright (C) 2013 Jonathan Currie (I2C2)

% This file will help you compile DSDP for use with MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get DSDP
% DSDP is available from http://www.mcs.anl.gov/hs/software/DSDP/. Download
% it and unzip the directory.

% 2) Compile DSDP
% The easiest way to compile DSDP is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer (you will need Intel MKL):
%
% %% Visual Studio Builder Commands
% path = 'full path to DSDP here'; %e.g. 'C:\Solvers\DSDP'
% sdir = [path '\src'];
% hdrs = [path '\include'];
% name = 'libdsdp';
% opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};
% VS_WriteProj(sdir,name,hdrs,opts)
% %%
%
% Once complete, you will have a directory called DSPD5.8\libdsdp. Open the
% Visual Studio solution, then complete the following steps:
%   a) Build a Win32 or x64 Release to compile the code.
%   b) Copy the generated .lib file to the following folder:
%
%   OPTI/Solvers/dsdp/Source/lib/win32 or win64
%
%   And rename it to libdsdp.lib. 
%
%   You will also need to copy all header files from DSDP/include to 
%   the following folder:
%
%   OPTI/Solvers/dsdp/Source/Include

% 3) DSDP MEX Interface
% The DSDP MEX Interface is a simple MEX interface I wrote to use DSDP 
% based in parts on the interface supplied with DSDP.

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the DSDP MEX file. Once you have completed all the
% above steps, simply run this file to compile DSDP! You MUST BE in the 
% base directory of OPTI!

clear dsdp

% Modify below function if it cannot find Intel MKL on your system.
mkl_link = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('DSDP MEX FILE INSTALL\n\n');

%Get Libraries
post = [' -IInclude -L' libdir ' -llibdsdp -llibut -output dsdp'];
%Get MKL Libraries (for BLAS)
post = [post mkl_link];

%CD to Source Directory
cdir = cd;
cd 'Solvers/dsdp/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims dsdpmex.c';
try
    eval([pre post])
    movefile(['dsdp.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:dsdp','Error Compiling DSDP!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');