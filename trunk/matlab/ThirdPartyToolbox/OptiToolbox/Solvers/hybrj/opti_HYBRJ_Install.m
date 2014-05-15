%% HYBRJ + HYBRD Install for OPTI Toolbox
% Copyright (C) 2012 Jonathan Currie (I2C2)

% This file will help you compile MINPACK HYBRJ + HYBRJ for use with MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get HYBRJ + HYBRD
% HYBRJ + HYBRD are part of MINPACK available from 
% http://www.netlib.org/minpack/. Download hybrd.f plus dependences and
% hybrj.f and unzip them all to a common directory. This way we will build
% one library with both hybrd and hybrj routines.

% 2) Compile HYBRJ + HYBRD
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
%   OPTI/Solvers/hybrj/Source/lib/win32 or win64
%
%   And rename it to libhybrj.lib. 

% 3) HYBRJ MEX Interface
% The HYBRJ MEX Interface is a simple MEX interface I wrote to use HYBRJ.

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the HYBRJ MEX file. Once you have completed all the
% above steps, simply run this file to compile HYBRJ! You MUST BE in the 
% base directory of OPTI!

clear hybrj

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('HYBRJ + HYBRD MEX FILE INSTALL\n\n');

%Get HYBRJ + HYBRD Libraries
post = [' -IInclude -L' libdir ' -llibhybrj -llibut -output hybrj'];
%Get Intel Fortran Libraries (for HYBRJ build) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link];

%CD to Source Directory
cdir = cd;
cd 'Solvers/hybrj/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims hybrjmex.c';
try
    eval([pre post])
    movefile(['hybrj.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:lmder','Error Compiling HYBRJ / HYBRD!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');