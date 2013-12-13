%% GLPK Install for OPTI Toolbox
% Copyright (C) 2012 Jonathan Currie (I2C2)

% This file will help you compile GNU Linear Programming Kit (GLPK) for use 
% with MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012

% To recompile you will need to get / do the following:

% 1) Get GLPK
% GLPK is available from http://www.gnu.org/software/glpk/glpk.html. We 
% will create VS projects below.

% 2) Compile GLPK
% The easiest way to compile GLPK is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:
%
% %% Visual Studio Builder Commands
% path = 'full path to GLPK here'; %e.g. 'C:\Solvers\GLPK'
% sdir = [path '\src'];
% name = 'libglpk';
% opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};
% VS_WriteProj(sdir,name,[],opts)
% %%
%
% Once complete, you will have a directory called GLPK/libglpk. Open the
% Visual Studio 2012 project file, then complete the following steps:
%   a) Build a Win32 or x64 Release to compile the code.
%   b) Copy the generated .lib file to the following folder:
%
%   OPTI/Solvers/glpk/Source/lib/win32 or win64
%
%   You will also need to copy all header files from GLPK/src to 
%   the following folder:
%
%   OPTI/Solvers/glpk/Source/Include

% 3) Get GLPKMEX
% GLPKMEX is a MEX interface to GLPK by Nicolo Giorgetti, and is used
% within OPTI to interface to GLPK. It is available from
% http://glpkmex.sourceforge.net/. Unzip the source files and copy
% glpkcc.cpp to:
%
%   OPTI/Solvers/glpk/Source

% Note I have made the following changes to the GLPKMEX file in the 
% supplied version:
% - Added Ctrl-C detection within print handler
% - Added drawnow to enable iteration by iteration printing
% - Changed the internal function description to be called 'glpk'

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the GLPK MEX file. Once you have completed all the
% above steps, simply run this file to compile GLPK! You MUST BE in the 
% base directory of OPTI!

clear glpk

% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('GLPK MEX FILE INSTALL\n\n');

%Get GLPK Libraries
post = [' -IInclude -L' libdir ' -llibglpk -llibut -output glpk'];

%CD to Source Directory
cdir = cd;
cd 'Solvers/glpk/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims glpkcc.cpp';
try
    eval([pre post])
    movefile(['glpk.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:glpk','Error Compiling GLPK!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');



