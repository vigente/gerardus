%% LP_SOLVE Install for OPTI Toolbox
% Copyright (C) 2012 Jonathan Currie (I2C2)

% This file will help you compile LP_SOLVE for use with MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012

% To recompile you will need to get / do the following:

% 1) Get LP_SOLVE
% LP_SOLVE is available from
% http://sourceforge.net/projects/lpsolve/files/lpsolve/. Download the
% source (lp_solve_5.5.2.0_source.tar.gz) or later version.

% 2) Compile LP_SOLVE
% The easiest way to compile LP_SOLVE is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the 
% required path on your computer:
%
% %% Visual Studio Builder Commands
% path = 'full path to LP_SOLVE here'; %e.g. 'C:\Solvers\LP_SOLVE'
% sdir = path;
% hdrs = [path '\bfp\bfp_LUSOL\LUSOL'];
% name = 'liblpsolve';
% opts.exPP = {'PARSER_LP','INVERSE_ACTIVE=INVERSE_LUSOL','RoleIsExternalInvEngine','_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE'};
% opts.exclude = {'lp_solveDLL.c'};
% opts.exFilter = {'lp_BFP*'};
% opts.exFolder = {'lpsolve55','lp_solve','demo','bfp\bfp_LUSOL\LUSOL'};
% VS_WriteProj(sdir,name,hdrs,opts)
% %%
%
% Once complete, you will have a directory called LP_SOLVE\liblpsolve. Open the
% Visual Studio 2012 project file, then complete the following steps:
%   a) Right click the project in the solution explorer, and click
%   Add -> Existing Item. Navigate to LP_SOLVE/bfp/bfp_LUSOL/LUSOL and add
%   lusol.c.
%   b) Build a Win32 or x64 Release to compile the code.
%   c) Copy the generated .lib file to the following folder:
%
%   OPTI/Solvers/lp_solve/Source/lib/win32 or win64
%
%   You will also need to copy all header files from LP_SOLVE to 
%   the following folder:
%
%   OPTI/Solvers/lp_solve/Source/Include

% 3) Get LP_SOLVE MEX Interface
% The LP_SOLVE MEX interface is available from the above link and is
% written by Peter Notebaert. Download the source
% (lp_solve_5.5.2.0_MATLAB_source.tar.gz) or later version. Copy the
% source files to the following folder:
%
%   OPTI/Solvers/lp_solve/Source
%
%   - lpsolve.c
%   - matlab.c

% Copy the header files to the following folder:
%
%   OPTI/Solvers/lp_solve/Source/Include
%
%   - lpsolvecaller.h
%   - matlab.h
%
% Note you DO NOT NEED the hash files or lp_explicit header.

% Make one change to matlab.h on line 12 to:
%
%   #if 1 (original = #if 0)
%
% Which will enable static linking (appears disabled by default)

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the LP_SOLVE MEX file. Once you have completed all 
% the above steps, simply run this file to compile LP_SOLVE! You MUST BE in 
% the base directory of OPTI!

clear lp_solve

% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('LP_SOLVE MEX FILE INSTALL\n\n');

%Get LP_SOLVE Libraries
post = [' -IInclude -L' libdir ' -lliblpsolve -DMATLAB -DWIN32 -DLPSOLVEAPIFROMLIB -output lp_solve'];

%CD to Source Directory
cdir = cd;
cd 'Solvers/lp_solve/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims lpsolve.c matlab.c';
try
    eval([pre post])
    movefile(['lp_solve.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:lpsolve','Error Compiling LP_SOLVE!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');