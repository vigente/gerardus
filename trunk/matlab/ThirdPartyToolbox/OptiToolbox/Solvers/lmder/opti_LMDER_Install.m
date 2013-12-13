%% LM_DER + LM_DIF Install for OPTI Toolbox
% Copyright (C) 2012 Jonathan Currie (I2C2)

% This file will help you compile MINPACK LM_DER + LM_DIF for use with MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get LM_DER + LM_DIF
% LM_DER + LM_DIF are part of MINPACK available from 
% http://www.netlib.org/minpack/. Download lmdif.f plus dependences and
% lmder.f and unzip them all to a common directory. This way we will build
% one library with both lm_dif and lm_der routines.

% 2) Compile LM_DER + LM_DIF
% The downloads contain multiple FORTRAN files which we will compile below.
% Follow the steps below to compile the library:
%   a) Create a new VS2012 Visual Fortran Static Library project.
%   b) Copy the FORTRAN files you downloaded to the project directory.
%   c) From the solution explorer right click Source Files and goto Add ->
%   Exisiting Item, and add all files.
%   d) Right click the project in the solution explorer, and click
%   properties. Navigate to Configuration Properties -> Fortran ->
%   Libraries and under Runtime Library change it to "Multithreaded DLL".
%   e) Create a new solution platform for x64 if required, ensure you are
%   building a 'Release', and build the project!
%   f) Copy the generated .lib file to the following folder:
%
%   OPTI/Solvers/lmder/Source/lib/win32 or win64
%
%   And rename it to liblmder.lib. 

% 3) LMDER MEX Interface
% The LMDER MEX Interface is a simple MEX interface I wrote to use LMDER.

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the LMDER MEX file. Once you have completed all the
% above steps, simply run this file to compile LMDER! You MUST BE in the 
% base directory of OPTI!

clear lmder

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('LM_DER + LM_DIF MEX FILE INSTALL\n\n');

%Get LMDER + LMDIF Libraries
post = [' -IInclude -L' libdir ' -lliblmder -llibut -output lmder'];
%Get Intel Fortran Libraries (for LMDER build) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link];

%CD to Source Directory
cdir = cd;
cd 'Solvers/lmder/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims lmdermex.c';
try
    eval([pre post])
    movefile(['lmder.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:lmder','Error Compiling LM_DER / LM_DIF!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');