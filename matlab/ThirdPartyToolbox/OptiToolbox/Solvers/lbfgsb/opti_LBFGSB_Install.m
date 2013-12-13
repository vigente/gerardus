%% L-BFGS-B Install for OPTI Toolbox
% Copyright (C) 2012 Jonathan Currie (I2C2)

% This file will help you compile Limited Memory Broyden-Fletcher-Goldfarb-
% Shanno Bounded Optimization (L-BFGS-B) for use with MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get L-BFGS-B
% L-BFGS-B is available from 
% http://users.eecs.northwestern.edu/~nocedal/lbfgsb.html. We will create 
% VS projects below.

% 2) Compile L-BFGS-B
% Version 3.0 contains multiple FORTRAN files which we will compile below.
% Follow the steps below to compile the library:
%   a) Create a new VS2012 Visual Fortran Static Library project.
%   b) Copy the FORTRAN files you downloaded to the
%   project directory.
%   c) From the solution explorer right click Source Files and goto Add ->
%   Exisiting Item, and add the following files:
%       - lbfgsb.f
%       - linpack.f
%       - timer.f
%       - blas.f (only if not compiling with MKL below)
%   d) Right click the project in the solution explorer, and click
%   properties. Navigate to Configuration Properties -> Fortran ->
%   Libraries and under Runtime Library change it to "Multithreaded DLL".
%   e) Create a new solution platform for x64 if required, ensure you are
%   building a 'Release', and build the project!
%   f) Copy the generated .lib file to the following folder:
%
%   OPTI/Solvers/lbfgsb/Source/lib/win32 or win64
%
%   And rename it to libLBFGSB.lib. 

% 3) Get L-BFGS-B MEX Interface
% Peter Carbonetto has written a MEX Interface to this solver and is
% available from http://www.cs.ubc.ca/~pcarbo/lbfgsb-for-matlab.html. I
% have modified this version slightly to incorporate OPTI functionality
% (function handles, easier options, and code changes for VC++ 2012), so I
% suggest you use the code supplied with this toolbox!

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the L-BFGS-B MEX file. Once you have completed all the
% above steps, simply run this file to compile L-BFGS-B! You MUST BE in the 
% base directory of OPTI!

clear lbfgsb

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('L-BFGS-B MEX FILE INSTALL\n\n');

%Get LBFGSB Libraries
post = [' -IInclude -L' libdir ' -llibLBFGSB -llibut'];
%Get Intel Fortran Libraries (for LBFGSB build) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link];

%CD to Source Directory
cdir = cd;
cd 'Solvers/lbfgsb/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims lbfgsb.cpp lbfgsb_program.cpp program.cpp';
try
    eval([pre post])
    movefile(['lbfgsb.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:lbfgsb','Error Compiling L-BFGS-B!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');