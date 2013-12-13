%% CSDP Install for OPTI Toolbox
% Copyright (C) 2013 Jonathan Currie (I2C2)

% This file will help you compile CSDP for use with MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get CSDP
% CSDP is available from https://projects.coin-or.org/Csdp/. Download
% it and unzip the directory.

% 2) Compile CSDP
% The easiest way to compile CSDP is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

% %% Visual Studio Builder Commands
% path = 'full path to CSDP here'; %e.g. 'C:\Solvers\CSDP'
% sdir = [path '\lib'];
% hdrs = [path '\include'];
% name = 'libcsdp';
% opts.exPP = {'_CRT_SECURE_NO_WARNINGS','USEOPENMP','NOSHORTS'};
% opts.openMP = true;
% opts.exclude = 'user_exit.c';
% VS_WriteProj(sdir,name,hdrs,opts)
% %%

% Once complete, you will have a directory called CSDP\libcsdp. Open the
% Visual Studio 2012 project file, then complete the following steps:
%   a) In order to accommodate more features of CSDP from MEX I modified 
%   sdp.c and easysdp.c. I suggest you copy my files from CSDP_mex_changes.zip
%   and replace the respective files in the CSDP lib directory. Also replace
%   declarations.h in the include directory. Changes include (for interest):
%       - easysdp.c:
%           - Adding struct paramstruc params as input arg
%           - Set printlevel = 0;
%           - Adding ppinf,pdinf,prealgap,pxzgap as extra return args
%           - Modified code to accept, return above appropriately
%       - sdp.c:
%           - Changed user_exit call to accept any non zero return arg
%           - Replaced all exit() calls with return instead
%           - Added two extra return codes for C and constraint checks
%       - declarations.h
%           - Updated function prototypes as per the above
%           - Declared malloc to _aligned_malloc so we can align memory for
%           SSE2 and BLAS/LAPACK (16 byte alignment)
%   b) Build a Win32 or x64 Release to compile the code.
%   c) Copy the generated .lib file to the following folder:
%
%   OPTI/Solvers/csdp/Source/lib/win32 or win64
%
%   You will also need to copy all header files from CSDP/include to 
%   the following folder:
%
%   OPTI/Solvers/csdp/Source/Include

% 3) CSDP MEX Interface
% The CSDP MEX Interface is a simple MEX interface I wrote to use CSDP.

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the CSDP MEX file. Once you have completed all the
% above steps, simply run this file to compile CSDP! You MUST BE in the 
% base directory of OPTI!

clear csdp

% Modify below function if it cannot find Intel MKL on your system.
mkl_link = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('CSDP MEX FILE INSTALL\n\n');

%Get Libraries
post = [' -IInclude -L' libdir ' -llibcsdp -llibut -output csdp'];
%Get MKL Libraries (for BLAS)
post = [post mkl_link];


%CD to Source Directory
cdir = cd;
cd 'Solvers/csdp/Source';

%Compile & Move (NOTE don't link against default VC++ OpenMP Lib as requires DLLs)
pre = 'mex -v -largeArrayDims LINKFLAGS="$LINKFLAGS /NODEFAULTLIB:vcomp.lib" -DNOSHORTS csdpmex.c';
try
    eval([pre post])
    movefile(['csdp.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:csdp','Error Compiling CSDP!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');