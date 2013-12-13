%% NL2SOL Install for OPTI Toolbox
% Copyright (C) 2012 Jonathan Currie (I2C2)

% This file will help you compile NL2SOL for use with MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get NL2SOL
% NL2SOL is available in multiple variants, however the most recent is 
% included in the PORT library, available from:
% http://netlib.sandia.gov/cgi-bin/netlib/netlibfiles.tar?filename=netlib/port
% The library contains functions for a range of mathematical functions,
% however we will just be using the NL2SOL variants.

% 2) Compile NL2SOL + NL2SNO (DN2F, DN2G, DN2FB, DN2GB)
% Follow the steps below to compile the library:
%   a) Create a new VS2010 Visual Fortran Static Library project.
%   b) From the solution explorer right click Source Files and goto Add ->
%   Exisiting Item, and add all .f files listed in Makefile. Note d1mach,
%   r1mach and i1mach.f should be added from the Mach/ folder.
%   c) We need to set the machine constants in function d1mach.f (and the
%   others if you intend on compiling single precision variants). The
%   default is IEEE Big Endian, which is not common these days. Comment
%   these out, and uncomment constants under IEEE 8087 (LSB - Little
%   Endian).
%   d) Right click the project in the solution explorer, and click
%   properties. Navigate to Configuration Properties -> Fortran ->
%   Libraries and under Runtime Library change it to "Multithreaded DLL".
%   e) Create a new solution platform for x64 if required, ensure you are
%   building a 'Release', and build the project!
%   f) Copy the generated .lib file to the following folder:
%
%   OPTI/Solvers/nl2sol/Source/lib/win32 or win64
%
%   And rename it to libnl2sol.lib. 

% 3) NL2SOL MEX Interface
% The NL2SOL MEX Interface is a simple MEX interface I wrote to use NL2SOL.

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the NL2SOL MEX file. Once you have completed all the
% above steps, simply run this file to compile NL2SOL! You MUST BE in the 
% base directory of OPTI!

clear nl2sol

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('NL2SOL MEX FILE INSTALL\n\n');

%Get NL2SOL Libraries
post = [' -IInclude -L' libdir ' -llibnl2sol -llibut -output nl2sol'];
%Get Intel Fortran Libraries (for NL2SOL build) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link];

%CD to Source Directory
cdir = cd;
cd 'Solvers/nl2sol/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims nl2solmex.c';
try
    eval([pre post])
    movefile(['nl2sol.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:nl2sol','Error Compiling NL2SOL!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');