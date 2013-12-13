%% IPOPT Install for OPTI Toolbox
% Copyright (C) 2012 Jonathan Currie (I2C2)

% This file will help you compile Interior Point OPTimizer (IPOPT) for use 
% with MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% NOTE - From OPTI v1.71 IPOPT is now dynamically linked against the
% MathWorks supplied libmwma57.dll (HSL MA57). This step is optional
% as we are also compiling MUMPS. Alternatively you can skip MUMPS, 
% and just use MA57! Also be aware libmwma57.dll does not play well
% on unconstrained problems due in part to missing MeTiS, thus the 
% ma57 pivot order option is overidden automatically.

% To recompile you will need to get / do the following:

% 1) Get IPOPT
% IPOPT is available from http://www.coin-or.org/download/source/Ipopt/.
% Download the latest version of the source.

% 2) Compile IPOPT
% The easiest way to compile IPOPT is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

% %% Visual Studio Builder Commands
% path = 'full path to IPOPT here'; %e.g. 'C:\Solvers\IPOPT'
% mumpspath = 'full path to MUMPS here'; %e.g. 'C:\Solvers\MUMPS'
% sdir = [path '\src'];
% hdrs = {[mumpspath '\include'], [mumpspath '\libseq'], [path '\..\BuildTools\headers']};
% name = 'libipopt';
% opts = [];
% opts.exPP = {'IPOPT_BUILD','FUNNY_MA57_FINT','COINHSL_HAS_MA57','_CRT_SECURE_NO_WARNINGS'};
% opts.exFolder = {'Apps','Algorithm\Inexact','contrib\LinearSolverLoader'};
% opts.exclude = 'AmplTNLP.cpp';
% VS_WriteProj(sdir,name,hdrs,opts)
% %%

% Once complete, you will have a directory called IPOPT\libipopt. Open the
% Visual Studio 2012 project file, then complete the following steps:
%   a) Open config_default.h and make the following changes:
%       - Add #define COIN_HAS_METIS 1 under the METIS library line
%       - Add #define COIN_HAS_MUMPS 1 under the MUMPS library line
%       - Comment the line #define COIN_HAS_HSL 1
%       - Comment the line #define COIN_HAS_ASL 1
%   b) Build a Win32 or x64 Release to compile the code.
%   c) Copy the generated .lib file to the following folder:
%
%   OPTI/Solvers/ipopt/Source/lib/win32 or win64
%
%   You will also need to copy all header files from IPOPT/src to 
%   the following folder (inc sub-directories):
%
%   OPTI/Solvers/ipopt/Source/Include

% 3) Complete MUMPS Compilation
% Supplied with this toolbox is the MUMPS linear system solver and METIS
% libraries. These are required when linking IPOPT so make sure you have
% built these first! See opti_MUMPS_Install.m in Solvers/mumps for more
% information.

% 4) libmwma57 Compilation
% As MathWorks only supply a dll, I have manually created an import library
% for libmwma57.dll, so we can link against it. The project and source for
% this is included in libmwma57.zip. It is optional to recompile this, as I
% have included compiled libraries for you!

% 5) IPOPT MEX Interface
% Supplied with the IPOPT source is a MEX interface by Dr. Peter Carbonetto.
% This is in Ipopt/contrib/MatlabInterface/src. However I have modified
% parts of this interface, so you must use the version supplied with OPTI.

% 6) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the IPOPT MEX file. Once you have completed all the
% above steps, simply run this file to compile IPOPT! You MUST BE in the 
% base directory of OPTI!

clear ipopt

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();
% Dependency Paths
mumpsdir = '..\..\mumps\Source\';

fprintf('\n------------------------------------------------\n');
fprintf('IPOPT MEX FILE INSTALL\n\n');

%Get IPOPT Libraries
post = ' -IInclude -IInclude/Common -IInclude/BuildTools -IInclude/Interfaces -IInclude/LinAlg -IInclude/LinAlg/TMatrices -IInclude/Algorithm ';
post = [post ' -L' libdir ' -llibIpopt -llibmwma57 -DhaveMA57'];
%Get MUMPS Libraries
post = [post ' -I' mumpsdir '\Include -L' mumpsdir libdir '-llibdmumps_c -llibdmumps_f -llibseq_c -llibseq_f -llibmetis -llibpord'];
%Get Intel Fortran Libraries (for MUMPS build) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link];

%Common Args
post = [post ' -DIPOPT_BUILD -output ipopt'];   

%CD to Source Directory
cdir = cd;
cd 'Solvers/ipopt/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims matlabexception.cpp matlabfunctionhandle.cpp matlabjournal.cpp iterate.cpp ipoptoptions.cpp options.cpp sparsematrix.cpp callbackfunctions.cpp matlabinfo.cpp matlabprogram.cpp ipopt.cpp';
try
    eval([pre post])
    movefile(['ipopt.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:ipopt','Error Compiling IPOPT!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');