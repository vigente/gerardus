%% OOQP Install for OPTI Toolbox
% Copyright (C) 2012 Jonathan Currie (I2C2)

% This file will help you compile Objective Orientated Quadratic
% Programming (OOQP) for use with MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get OOQP
% OOQP is available from http://pages.cs.wisc.edu/~swright/ooqp/. You will 
% need to register before you can download. We will create the VS project 
% below.

% 2) Compile OOQP
% The easiest way to compile OOQP is to use the Visual Studio Project
% Builder included with OPTI. To begin with, copy the files in
% ooqp/Source/Pardiso Addin.zip to a suitable folder under the /src/
% directory in OOQP. This will add my PARDISO Add In (optional step).

%Use the following commands, substituting the required path on your computer:
%
% %% Visual Studio Builder Commands
% path = 'full path to OOQP here'; %e.g. 'C:\Solvers\OOQP'
% sdir = [path '\src'];
% name = 'libooqp';
% opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};
% opts.mkl = true;
% opts.exclude = {'QpGenSparseOblio.C','QpGenSparseSuperLu.C','OoqpPetscMonitor.C',...
%                 'QpGenSparseGondzioDriver.C','QpGenSparseMa57GondzioDriver.C','QpGenSparseMehrotraDriver.C',...
%                 'QpBoundDenseGondzioDriver.C','QpGenDenseGondzioDriver.C','SvmGondzioDriver.C','HuberGondzioDriver.C'};
% opts.exFilter = {'QpBoundPetsc*'};
% opts.exFolder = {'Ampl','CInterface','PetscLinearAlgebra','QpExample','Mex'};
% VS_WriteProj(sdir,name,[],opts)
% %%
%
% Once complete, you will have a directory called OOQP\libooqp. Open the
% Visual Studio 2012 project file, then complete the following steps:
%   a) Right click the project in the solution explorer, and click
%   properties. Navigate to C/C++ -> Code Generation and under
%   Enable C++ Exceptions change to "Yes with Extern C Functions (/EHs)".
%   b) Navigate to C/C++ -> Advanced and under Compile As change to
%   "Compile as C++ Code (/TP)". This is required as the file extension is
%   .c, but we are compiling C++ code.
%   c) If you are building a 64bit build AND using MATLAB's MA57 then add
%   the preprocessor 'X64MA57' under C/C++ -> Preprocessor Definitions.
%   d) Replace Ma57Solver.h and Ma57Solver.C with those in  
%   Solvers/ooqp/Source/Ma57Solver.zip.
%   e) Open hash.C and on lines 109, 140 and 163 change "hash" to "hashL".
%   f) Change MehrotraSolver.C line 33 to "extern double gmu;" 
%   g) If you are going to compile with MA27 and Intel Fortran you
%   will need to rename the MA27 function calls. Open Ma27Solver.h and on
%   lines 16, 18, 26, and 35 rename the ma27 function calls as below:
%       - ma27id_  becomes  MA27ID
%       - ma27ad_  becomes  MA27AD, and so forth for ma27bd_ and ma27cd_.
%   You will also need to rename them in Ma27Solver.c so you may like to do
%   a find and replace on each one.
%   h) Build a Win32 or x64 Release to compile the code.
%   i) Copy the generated .lib file to the following folder:
%
%   OPTI/Solvers/ooqp/Source/lib/win32 or win64
%
%   You will also need to copy all header files from OOQP/src to 
%   the following folder:
%
%   OPTI/Solvers/ooqp/Source/Include

% 4) MEX Interface
% The MEX interface is a modification of the one supplied with OOQP.

% 5) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the OOQP MEX file. Once you have completed all the
% above steps, simply run this file to compile OOQP! You MUST BE in the 
% base directory of OPTI!

clear ooqp

%Change if linking against MA27 as well
haveMA27 = false;

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('OOQP MEX FILE INSTALL\n\n');

%Get OOQP Libraries
post = [' -IInclude -L' libdir ' -llibooqp -llibut -output ooqp'];
%Get MA57 Library
post = [post ' -llibmwma57 -DHAVE_MA57'];
%Include optional MA27
if(haveMA27)
	post = [post ' -lma27 -DHAVE_MA27'];
end
%Get Intel Fortran Libraries (for MA27 build if required) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link ' -DHAVE_PARDISO'];

%CD to Source Directory
cdir = cd;
cd 'Solvers/ooqp/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims ooqp_mex.cpp';
try
    eval([pre post])
    movefile(['ooqp.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:ooqp','Error Compiling OOQP!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');
