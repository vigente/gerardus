%% MUMPS Install for OPTI Toolbox
% Copyright (C) 2011 Jonathan Currie (I2C2)

% This file will help you compile aMUltifrontal Massively Parallel sparse
% direct Solver (MUMPS) for use with MATLAB. 

% The supplied files and instructions are for compiling sequential double 
% precision MUMPS only.

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get MUMPS
% MUMPS is available from http://graal.ens-lyon.fr/MUMPS/. You will need to
% register before you can download. We will create VS projects below.

% 2) Get METIS
% METIS is available from http://glaros.dtc.umn.edu/gkhome/fsroot/sw/metis/OLD.
% Download version 4.0.3 (last compatible version with MUMPS).

% 3) Compile MUMPS and METIS
% The easiest way to compile MUMPS is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

% %% Visual Studio Builder Commands
% clear
% path = 'full path to MUMPS here'; %e.g. 'C:\Solvers\MUMPS'
% metispath = 'full path to METIS here'; e.g. 'C:\Solvers\METIS'
% n = 1;
% % DMUMPS_C
% sdir = [path '\src'];
% hdrs = {[path '\include'],[path '\PORD\include']};
% name = 'libdmumps_c';
% opts = [];
% opts.exPP = {'_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE','MUMPS_ARITH=MUMPS_ARITH_d','pord','metis'};
% VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% % DMUMPS_F
% sdir = [path '\src'];
% hdrs = {[path '\include'],[path '\libseq']};
% name = 'libdmumps_f';
% opts = [];
% opts.cpp = false;
% opts.exPP = {'pord','metis'};
% opts.exFilter = {'cmumps*','smumps*','zmumps*'};
% VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% % ZMUMPS_C 
% sdir = [path '\src'];
% hdrs = {[path '\include'],[path '\PORD\include']};
% name = 'libzmumps_c';
% opts = [];
% opts.exPP = {'_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE','MUMPS_ARITH=MUMPS_ARITH_z','pord','metis'};
% VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% % ZMUMPS_F 
% sdir = [path '\src'];
% hdrs = {[path '\include'],[path '\libseq']};
% name = 'libzmumps_f';
% opts = [];
% opts.cpp = false;
% opts.exPP = {'pord','metis'};
% opts.exFilter = {'cmumps*','smumps*','dmumps*'};
% VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% % LIBSEQ_C
% sdir = [path '\libseq'];
% name = 'libseq_c';
% opts = [];
% opts.exPP = {'_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE'};
% VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% % LIBSEQ_F
% sdir = [path '\libseq'];
% name = 'libseq_f';
% opts = [];       
% opts.cpp = false;
% VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% % PORD
% sdir = [path '\PORD\lib'];
% hdrs = [path '\PORD\include'];
% name = 'libpord';
% opts = [];
% opts.exPP = {'_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE'};
% VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% % METIS
% sdir = [metispath '\Lib'];
% name = 'libmetis';
% opts = [];
% opts.exPP = {'_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE','__STDC__','__VC__'};
% VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% %Write the Solution File
% VS_WriteSol(VSPRJ)
% %%

% Once complete, you will have a directory called MUMPS/libdmumps_c. Open 
% the Visual Studio 2012 solution file, then complete the following steps:
%   a) In the project 'libmetis' complete the following steps:
%   	- Open metis.h and comment out line 21 (#include strings.h)
%   	- Open proto.h and uncomment line 435 (void GKfree...)
%       - (METIS >= v5 only) in metis.h define REALTYPEWIDTH 64
%       - (METIS >= v5 only) in metis.h define USE_GKREGEX
%   b) Build a Win32 or x64 Release of each proejct to compile the code.
%   c) Copy the generated .lib files to the following folder:
%
%   OPTI/Solvers/mumps/Source/lib/win32 or win64
%
%   You will also need to copy all header files from MUMPS/src 
%   to the following folder:
%
%   OPTI/Solvers/mumps/Source/Include

% 4) Get MUMPS MEX Interface
% Supplied with MUMPS in the MATLAB folder is mumpsmex.c. Copy this to the
% following folder:
%
%   OPTI/Solvers/mumps/Source 

% 5) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the MUMPS MEX file. Once you have completed all the
% above steps, simply run this file to compile MUMPS! You MUST BE in the 
% base directory of OPTI!

clear mumps

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('MUMPS MEX FILE INSTALL\n\n');

%Get MUMPS Libraries
post = [' -IInclude -L' libdir ' -llibdmumps_c -llibdmumps_f -llibpord -llibseq_c -llibseq_f -llibmetis -DMUMPS_ARITH=2'];
%Get Intel Fortran Libraries (for MUMPS build) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link];
%Common
post = [post ' -output mumps'];

%CD to Source Directory
cdir = cd;
cd 'Solvers/mumps/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims mumpsmex.c';
try
    eval([pre post])
    movefile(['mumps.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:mumps','Error Compiling MUMPS!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');

% METIS Reference:
% “A Fast and Highly Quality Multilevel Scheme for Partitioning Irregular 
% Graphs”. George Karypis and Vipin Kumar. SIAM Journal on Scientific 
% Computing, Vol. 20, No. 1, pp. 359—392, 1999.

%% Optional ZMUMPS Install
clear zmumpsmex

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('ZMUMPS MEX FILE INSTALL\n\n');

%Get MUMPS Libraries
post = [' -IInclude -L' libdir ' -llibzmumps_c -llibzmumps_f -llibpord -llibseq_c -llibseq_f -llibmetis -DMUMPS_ARITH=8'];
%Get Intel Fortran Libraries (for MUMPS build) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link];
%Common
post = [post ' -output zmumpsmex'];

%CD to Source Directory
cdir = cd;
cd 'Solvers/mumps/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims mumpsmex.c';
try
    eval([pre post])
    movefile(['zmumpsmex.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:mumps','Error Compiling ZMUMPS!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');