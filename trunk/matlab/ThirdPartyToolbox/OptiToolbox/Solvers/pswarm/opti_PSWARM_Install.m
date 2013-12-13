%% PSWARM Install for OPTI Toolbox
% Copyright (C) 2012 Jonathan Currie (I2C2)

% This file will help you compile PSwarm for use with MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 1) Get PSwarm
% PSwarm is available from http://www.norg.uminho.pt/aivaz/pswarm/. 
% Download PPSwarm_vxx.zip (C Version) and unzip to a suitable location.

% 2) Compile PSwarm
% The easiest way to compile PSwarm is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the 
% required path on your computer:
%
% %% Visual Studio Builder Commands
% path = 'full path to PSwarm here'; %e.g. 'C:\Solvers\PSwarm'
% sdir = path;
% name = 'libpswarm';
% opts.exPP = {'LINEAR','_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE'};
% opts.exclude = {'showcache.c','user.c','pswarm_main.c','pswarm_py.c','pswarm_r.c','cache.c'};
% VS_WriteProj(sdir,name,[],opts)
% %%
%
% Once complete, you will have a directory called PSwarm\libpswarm. Open the
% Visual Studio 2012 project file, then complete the following steps:
%   a) A couple code changes are required in order to interface with OPTI:
%       - In pswarm.h add "int solveriters;" to the Stats structure as the
%       last line.
%       - In pswarm.h change the define on line 62 to 
%         "#define SYS_RANDOM 1"
%       - In pswarm.c on line 666 add "stats.solveriters = iter;"
%       - Change the structure definitions in pswarm.c at line 66 to extern:
%           "extern struct Stats stats;"
%           "extern struct Options opt;"
%       - Comment the default options structure in pswarm.c on line 70
%       (comment all lines for this structure).
%   b) Build a Win32 or x64 Release to compile the code.
%   c) Copy the generated .lib file to the following folder:
%
%   OPTI/Solvers/pswarm/Source/lib/win32 or win64
%
%   You will also need to copy pswarm.h to the following folder:
%
%   OPTI/Solvers/pswarm/Source/Include

% 3) PSwarm MEX Interface
% The PSwarm MEX Interface is a simple MEX interface I wrote to use PSwarm. 

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the PSwarm MEX file. Once you have completed all the
% above steps, simply run this file to compile PSwarm! You MUST BE in the 
% base directory of OPTI!

clear pswarm

% Modify below function if it cannot find Intel MKL on your system.
mkl_link = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('PSwarm MEX FILE INSTALL\n\n');

%Get PSwarm Libraries
post = [' -IInclude -L' libdir ' -llibpswarm -llibut -output pswarm'];
%Get MKL Libraries (for BLAS)
post = [post mkl_link];

%CD to Source Directory
cdir = cd;
cd 'Solvers/pswarm/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims pswarmmex.c';
try
    eval([pre post])
    movefile(['pswarm.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:pswarm','Error Compiling PSwarm!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');