%% FILTERSD Install for OPTI Toolbox
% Copyright (C) 2013 Jonathan Currie (I2C2)

% This file will help you compile FILTER SD for use with MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% 1) Get FILTER SD
% FilterSD is available from http://www.coin-or.org/download/source/filterSD/
% Download the .zip file for use on Windows.

% 2) Compile FILTERSD + FILTERSDSP
% The downloads contain multiple FORTRAN files which we will compile below.
% Follow the steps below to compile the library:
%   a) Create a new VS2010 Visual Fortran Static Library project.
%   b) Copy the FORTRAN files you downloaded to the project directory.
%   c) From the solution explorer right click Source Files and goto Add ->
%   Exisiting Item, and add the following files:
%       - checkd.f
%       - denseA.f
%       - denseL.f
%       - filterSD.f
%       - glcpd.f
%       - l1sold.f
%       - util.f
%
%   d) Replace filterSD.f and glcpd.f with the files located in the
%   following folder (only compatible with OPTI version):
%       OPTI/Solvers/filtersd/Source/filterSD_JCEdit.zip
%   (They contain an edit to allow extra solver termination conditions)  
%   e) Right click the project in the solution explorer, and click
%   properties. Navigate to Configuration Properties -> Fortran ->
%   Libraries and under Runtime Library change it to "Multithreaded DLL".
%   f) Create a new solution platform for x64 if required, ensure you are
%   building a 'Release', and build the project!
%   g) Copy the generated .lib file to the following folder:
%
%   OPTI/Solvers/filtersd/Source/lib/win32 or win60
%
%   And rename it to libfilterSD.lib. 
%
%   h) Under the same solution, follow steps a-g again, but this time add
%   the following files to the project:
%       - checkd.f
%       - filterSD.f
%       - glcpd.f
%       - l1sold.f
%       - schurQR.f
%       - sparseA.f
%       - util.f
%
%   and copy and rename the generated library to libfilterSDsp.lib.

% 3) FILTERSD MEX Interface
% The FILTERSD MEX Interface is a simple MEX interface I wrote to use
% FILTERSD.

% 6) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the FILTERSD MEX file. Once you have completed all the
% above steps, simply run this file to compile FILTERSD! You MUST BE in the 
% base directory of OPTI!

clear filtersd filtersdsp

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('FILTERSD MEX FILE INSTALL\n\n');

%Get FILTERSD Libraries
post = [' -L' libdir ' -llibfilterSD -llibut'];
%Get Intel Fortran Libraries (for ifort build) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link];
%Common Args
post = [post ' -output filtersd'];   

%CD to Source Directory
cdir = cd;
cd 'Solvers/filtersd/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims filtersdmex.c';
try
    eval([pre post])
    movefile(['filtersd.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:filtersd','Error Compiling FILTERSD!\n%s',ME.message);
end

%Compile & Move Sparse Version
post = regexprep(post,'-llibfilterSD','-llibfilterSDsp');
post = regexprep(post,'-output filtersd','-output filtersdsp');
try
    eval([pre post ' -DSPARSEVER'])
    movefile(['filtersdsp.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:filtersd','Error Compiling FILTERSD SPARSE!\n%s',ME.message);
end

cd(cdir);
fprintf('------------------------------------------------\n');