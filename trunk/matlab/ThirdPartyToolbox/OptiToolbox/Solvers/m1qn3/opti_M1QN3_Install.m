%% M1QN3 Install for OPTI Toolbox
% Copyright (C) 2012 Jonathan Currie (I2C2)

% This file will help you compile M1QN3 for use with MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get M1QN3
% M1QN3 is available from 
% https://who.rocq.inria.fr/Jean-Charles.Gilbert/modulopt/optimization-routines/m1qn3/m1qn3.html
% Download and Unzip the folder, we will create a Visual Studio project below.

% 2) Compile M1QN3
% Follow the steps below to compile the library:
%   a) Create a new VS2012 Visual Fortran Static Library project.
%   b) From the solution explorer right click Source Files and goto Add ->
%   Exisiting Item, and add m1qn3.f from src/.
%   c) Right click the project in the solution explorer, and click
%   properties. Navigate to Configuration Properties -> Fortran ->
%   Libraries and under Runtime Library change it to "Multithreaded DLL".
%   d) Create a new solution platform for x64 if required, ensure you are
%   building a 'Release', and build the project!
%   e) Copy the generated .lib file to the following folder:
%
%   OPTI/Solvers/m1qn3/Source/lib/win32 or win64
%
%   And rename it to libm1qn3.lib. 

% 3) M1QN3 MEX Interface
% The M1QN3 MEX Interface is a simple MEX interface I wrote to use M1QN3.

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the NL2SOL MEX file. Once you have completed all the
% above steps, simply run this file to compile M1QN3! You MUST BE in the 
% base directory of OPTI!

clear m1qn3

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('M1QN3 MEX FILE INSTALL\n\n');

%Get M1QN3 Libraries
post = [' -IInclude -L' libdir ' -llibm1qn3 -llibut -output m1qn3'];
%Get Intel Fortran Libraries (for M1QN3 build) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link];

%CD to Source Directory
cdir = cd;
cd 'Solvers/m1qn3/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims m1qn3mex.c';
try
    eval([pre post])
    movefile(['m1qn3.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:m1qn3','Error Compiling M1QN3!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');