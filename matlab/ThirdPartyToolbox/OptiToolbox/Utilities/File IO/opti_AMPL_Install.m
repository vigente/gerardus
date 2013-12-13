%% AMPL Install for OPTI Toolbox
% Supplied binaries are built from Netlib's AMPL Solver Library Interface

%   Copyright (C) 2011 Jonathan Currie (I2C2)

% This file will help you compile AMPL Solver Library (ASL) for use with 
% MATLAB.

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012

% 1) Get AMPL Solver Library
% The generic NL reader for AMPL is available free from Netlib 
% (http://www.netlib.org/ampl/solvers/). You will need to download all .c
% and .h as well as .hd files. Note this is not the AMPL engine
% (www.ampl.com) which is a commerical product, but code to allow people to
% connect their solvers to AMPL. Alternatively send a blank email to 
% "netlib@netlib.org" with "send all from ampl/solvers" as the subject to 
% retrieve all files.

% 2) Compile ASL
% The easiest way to compile ASL is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

% %% Visual Studio Builder Commands
% path = 'full path to ASL here'; %e.g. 'C:\Solvers\ASL'
% sdir = path;
% name = 'libasl';
% opts = [];
% opts.exPP = {'Arith_Kind_ASL=1','Sscanf=sscanf','Printf=printf','Sprintf=sprintf',...
%              'Fprintf=fprintf','snprintf=_snprintf','NO_STDIO1','_CRT_SECURE_NO_WARNINGS','_CRT_NONSTDC_NO_DEPRECATE'};
% opts.exclude = {'arithchk.c','atof.c','b_search.c','dtoa.c','fpinit.c','funcadd.c',...
%                 'funcadd0.c','funcaddk.c','funcaddr.c','obj_adj0.c','sjac0dim.c',...
%                 'sprintf.c','sscanf.c','printf.c','mpec_adj0.c'};
% opts.charset = 'MultiByte';            
% VS_WriteProj(sdir,name,[],opts)
% %%

% Once complete, you will have a directory called ASL\libasl. Open the
% Visual Studio 2012 project file, then complete the following steps:
%   a) Rename arith.h0 to arith.h, and uncomment the define "IEEE_8087"
%   b) Rename stdio1.h0 to stdio1.h
%   c) Under project properties change the character set to "Use
%   Multi-Byte Character Set" (both Win32 and x64 as required).
%   d) Comment out the line "exit(n)" in mainexit.c (this stops MATLAB 
%   crashing on an ASL error) [line 62]
%   e) To prevent Matlab crashing on a bad NL file comment all code under
%   the if(n_con < 0 || ...) statement on line 260 in jac0dim.c. Add return 
%   NULL; so instead of exiting, our MEX file can determine this as a bad file.
%   f) Comment line 900 in asl.h (extern int Sprintf(char*, ....)
%   g) Comment line 36 in jac0dim.c (extern int Sscanf(char*, ...)
%   h) Comment line 56 in stderr.c (AllocConsole();)
%   i) Build a Win32 or x64 Release to compile the code.
%   j) Copy the generated .lib file to the following folder:
%
%   OPTI/Utilities/File IO/Source/lib/win32 or win64
%
%   Also copy all header files to the include directory as per the above
%   path.

% 3) MEX Interface
% The Read MEX Interface is a simple MEX interface I wrote to use the AMPL 
% File IO and Eval routines.

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the MEX file. Once you have completed all 
% the above steps, simply run this file to compile! You MUST BE in 
% the base directory of OPTI!

clear asl

% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('AMPL MEX FILE INSTALL\n\n');

post = [' -IInclude -L' libdir ' -llibasl -DNO_STDIO1 -output asl'];

%CD to Source Directory
cdir = cd;
cd 'Utilities/File IO/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims amplmex.c';
try
    eval([pre post])
     movefile(['asl.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:ampl','Error Compiling AMPL!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');
