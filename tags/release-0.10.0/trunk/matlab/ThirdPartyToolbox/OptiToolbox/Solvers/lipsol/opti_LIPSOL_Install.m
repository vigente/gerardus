%% LIPSOL Install for OPTI Toolbox
% Copyright (C) 2013 Jonathan Currie (I2C2)

% This file will help you compile LIPSOL for use with MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012
% - Intel Compiler XE (FORTRAN)

% To recompile you will need to get / do the following:

% 1) Get LIPSOL
% LIPSOL is available from http://www.caam.rice.edu/~zhang/lipsol/.
% Download the latest version.

% 3) Compile LIPSOL FORTRAN Routines
% The easiest way to compile LIPSOL is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

% %% Visual Studio Builder Commands
% clear
% path = 'full path to LIPSOL here'; %e.g. 'C:\Solvers\LIPSOL'
% sdir = [path '\src'];
% name = 'liblipsol';
% opts = [];
% opts.cpp = false;
% opts.exclude = {'blkfctg.f','blkslvg.f','inpnvg.f','mps2mat.f','mps2matg.f','ordmmdg.f','symfctg.f'};
% VS_WriteProj(sdir,name,[],opts)
% %%

% Once complete, you will have a directory called LIPSOL/liblipsol. Open 
% the Visual Studio 2012 solution file, then complete the following steps:
%   a) Build a Win32 or x64 Release of each proejct to compile the code.
%   b) Copy the generated .lib files to the following folder:
%
%   OPTI/Solvers/lipsol/Source/lib/win32 or win64

% 3) LIPSOL MEX Interface
% LIPSOL can with a FORTRAN MEX Interface, however I have re-written it in
% C, adding 64bit Compatibility.

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the LIPSOL MEX files. Once you have completed all the
% above steps, simply run this file to compile LIPSOL! You MUST BE in the 
% base directory of OPTI!

clear ls_blkfct ls_blkslv ls_inpnv ls_ordmmd ls_symfct

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('LIPSOL FILE INSTALL\n\n');

%Files to Compile
files = {'ordmmdg','blkslvg','blkfctg','symfctg','inpnvg'};

%Get LIPSOL Libraries
post = [' -L' libdir ' -lliblipsol -llibut'];
%Get Intel Fortran Libraries (for LIPSOL build) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link];

%CD to Source Directory
cdir = cd;
cd 'Solvers/lipsol/Source';

%Compile each file & Move
pre = 'mex -v -largeArrayDims ';

for i = 1:length(files)
    try
        eval([pre files{i} '.c' post ' -output ls_' files{i}(1:end-1)])
        movefile(['ls_' files{i}(1:end-1) '.' mexext],'../','f')
        fprintf('Done!\n');
    catch ME
        cd(cdir);
        error('opti:lipsol','Error Compiling LIPSOL: %s!\n%s',files{i},ME.message);
    end
end
cd(cdir);
fprintf('------------------------------------------------\n');