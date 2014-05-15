%% NLOPT Install for OPTI Toolbox
% Copyright (C) 2012 Jonathan Currie (I2C2)

% This file will help you compile NonLinear OPTimization (NLOPT) for use 
% with MATLAB. 

% My build platform:
% - Windows 8 SP1 x64
% - Visual Studio 2012

% To recompile you will need to get / do the following:

% 1) Get NLOPT
% NLOPT is available from http://ab-initio.mit.edu/wiki/index.php/NLopt. 
% Download the source.

% 2) Compile NLOPT
% The easiest way to compile NLOPT is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the 
% required path on your computer:
%
% %% Visual Studio Builder Commands
% path = 'full path to NLOPT here'; %e.g. 'C:\Solvers\NLOPT'
% sdir = path;
% name = 'libnlopt';
% opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};
% opts.exclude = {'testfuncs.c','tst.cc','tstc.c','testros.cc','prog.cc','redblack_test.c','DIRparallel.c'};
% opts.exFolder = {'octave','swig','test'};
% VS_WriteProj(sdir,name,[],opts)
% %%
%
% Once complete, you will have a directory called NLOPT\libnlopt. Open the
% Visual Studio 2012 project file, then complete the following steps:
%   a) A default config.h has not been supplied (and I'm guessing only gets
%   made in Linux versions) so I've made one up based on the template. Copy
%   the supplied config.h from:
%       OPTI/Solvers/nlopt/Source/Include
%   To the project directory.
%   b) Build a Win32 or x64 Release to compile the code.
%   c) Copy the generated .lib file to the following folder:
%
%   OPTI/Solvers/nlopt/Source/lib/win32 or win64
%
%   You will also need to copy nlopt.h to the following folder:
%
%   OPTI/Solvers/nlopt/Source/Include
% 
%   Note with Global Optimization turned on in VS it will take a long time
%   to link the MEX file.

% 3) NLOPT MEX Interface
% The NLOPT MEX Interface was written by Steven Johnson and is located in
% the octave folder (nlopt_optimize-mex.c) HOWEVER in its original form it
% was not compatible with OPTI (due to the method of adding nonlinear
% constraints cell by cell), as well as the file name caused compile
% problems (no '-' allowed for MEX). Therefore I have modified the
% MEX interface and included is nlopt_optimize_mex.c which you will need to
% use for compatibility with OPTI. You will need to however copy
% nlopt_optimize_usage.h from nlopt-xxx/octave to:
%
%   OPTI/Solvers/nlopt/Source/Include

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the NLOPT MEX file. Once you have completed all 
% the above steps, simply run this file to compile NLOPT! You MUST BE in 
% the base directory of OPTI!

clear nlopt

% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('NLOPT MEX FILE INSTALL\n\n');

%Get NLOPT Libraries
post = [' -IInclude -L' libdir ' -llibnlopt -llibut -output nlopt'];

%CD to Source Directory
cdir = cd;
cd 'Solvers/nlopt/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims nloptmex.c';
try
    eval([pre post])
    movefile(['nlopt.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:nlopt','Error Compiling NLOPT!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');
