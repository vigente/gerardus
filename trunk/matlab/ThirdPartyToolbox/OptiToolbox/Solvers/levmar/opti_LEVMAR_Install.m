%% LEVMAR Install for OPTI Toolbox
% Copyright (C) 2012 Jonathan Currie (I2C2)

% This file will help you compile Levenberg-Marquardt in C/C++ (LEVMAR) for 
% use with MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get LEVMAR
% LEVMAR is available from http://www.ics.forth.gr/~lourakis/levmar/. We 
% will create the VS project below.

% 2) Compile LEVMAR
% The easiest way to compile LEVMAR is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the 
% required path on your computer:
%
% %% Visual Studio Builder Commands
% path = %'full path to LEVMAR here'; %e.g. 'C:\Solvers\LEVMAR'
% sdir = path;
% name = 'liblevmar';
% opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};
% opts.exclude = {'expfit.c','Axb_core.c','lm_core.c','lmbc_core.c','lmblec_core.c','lmbleic_core.c','lmlec_core.c','misc_core.c'};
% opts.exFolder = {'matlab'};
% VS_WriteProj(sdir,name,[],opts)
% %%
%
% Once complete, you will have a directory called LEVMAR\liblevmar. Open the
% Visual Studio 2012 project file, then complete the following steps:
%   a) Build a Win32 or x64 Release to compile the code.
%   b) Copy the generated .lib file to the following folder:
%
%   OPTI/Solvers/levmar/Source/lib/win32 or win64
%
%   You will also need to copy levmar.h to the following folder:
%
%   OPTI/Solvers/levmar/Source/Include

% 3) MEX Interface
% I have heavily modified the levmar supplied MEX interface thus I suggest
% you use my version, levmarmex.c, included as part of this distribution.
% This will also maintain compatibility with OPTI.

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the LEVMAR MEX file. Once you have completed all the
% above steps, simply run this file to compile LEVMAR! You MUST BE in the 
% base directory of OPTI!

clear levmar

% Modify below function if it cannot find Intel MKL on your system.
mkl_link = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('LEVMAR MEX FILE INSTALL\n\n');

%Get LEVMAR Libraries
post = [' -IInclude -L' libdir ' -lliblevmar'];
%Get MKL Libraries (for BLAS & LAPACK)
post = [post mkl_link];
%Common outputs
post = [post ' -output levmar'];

%CD to Source Directory
cdir = cd;
cd 'Solvers/levmar/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims levmarmex.c';
try
    eval([pre post])
    movefile(['levmar.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:levmar','Error Compiling LEVMAR!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');
