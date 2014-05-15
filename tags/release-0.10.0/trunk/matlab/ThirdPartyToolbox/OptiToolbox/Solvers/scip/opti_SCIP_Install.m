%% SCIP Install for OPTI Toolbox
% Copyright (C) 2012 Jonathan Currie (I2C2)

% This file will help you compile Solving Constraint Integer Programs (SCIP) 
% for use with MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 0) Complete Compilation as per OPTI instructions for MUMPS and IPOPT, in 
% that order.

% 1) Get SCIP & SoPlex
% SCIP is available from http://scip.zib.de/. Ensure you download the SCIP
% optimization suite, which includes SoPlex. We will create the VS projects 
% below.

% 2) Compile SCIP & SoPlex
% The easiest way to compile SCIP is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

% %% Visual Studio Builder Commands
% clear
% path = 'full path to SCIP here'; %e.g. 'C:\Solvers\SCIP'
% splxpath = 'full path to SOPLEX here'; %e.g. 'C:\Solvers\SOPLEX'
% ipoptpath = 'full path to IPOPT here'; %e.g. 'C:\Solvers\IPOPT'
% cppadpath = 'full path to CppAD here'; %e.g. 'C:\Solvers\CppAD'
% n = 1;
% % SCIP
% sdir = [path '\src'];
% hdrs = {[ipoptpath '\src'], [ipoptpath '\..\BuildTools\headers'], [splxpath '\src'], cppadpath};
% name = 'libscip';
% opts = [];
% opts.exPP = {'IPOPT_BUILD','_CRT_SECURE_NO_WARNINGS','NO_RAND_R','NO_SIGACTION','NO_STRERROR_R',...
%              'NO_STRTOK_R','NO_NEXTAFTER','ROUNDING_MS','NPARASCIP'};          
% opts.exclude = {'exprinterpret_none.c','nlpi_ipopt_dummy.c','nlpi_xyz.c','lpi_none.c',...
%                 'lpi_clp.cpp','lpi_cpx.c','lpi_grb.c','lpi_msk.c','lpi_qso.c',...
%                 'lpi_spx121.cpp','lpi_spx132.cpp','lpi_xprs.c','sorttpl.c','cmain.c',...
%                 'cppmain.cpp','disp_xyz.c','branch_xyz.c','event_xyz.c','cons_xyz.c',...
%                 'heur_xyz.c','presol_xyz.c','prop_xyz.c','pricer_xyz.c','nodesel_xyz.c',...
%                 'relax_xyz.c','reader_xyz.c','sepa_xyz.c','dialog_xyz.c'};
% VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% % SOPLEX
% sdir = [splxpath '\src'];
% name = 'libsoplex';
% opts = [];
% opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};
% opts.exclude = {'soplexmain.cpp','simpleexample.cpp'};
% VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = []; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% %Write the Solution File
% VS_WriteSol(VSPRJ)
% %%

% Once complete, you will have a directory called SCIP/libscip. Open the
% Visual Studio 2012 solution file, then complete the following steps:
%   a) (SCIP) Versions including and prior to 3.0.1 require line 2905 in
%   cons_soc.c to be changed to "if( lhscount >= nquadvars - 1 )"
%   b) (SCIP) Due to the way we have created the project the header files 
%   cannot be found. Perform a find and replace on the CURRENT PROJECT of the
%   following terms:
%       - FIND #include "tclique/       REPLACE #include "../tclique/
%       - FIND #include "blockmemshell/ REPLACE #include "../blockmemshell/
%       - FIND #include <blockmemshell/ REPLACE #include <../blockmemshell/
%       - FIND #include "scip/          REPLACE #include "../scip/
%       - FIND #include "nlpi/          REPLACE #include "../nlpi/
%       - FIND #include "objscip/       REPLACE #include "../objscip/
%       - FIND #include "xml/           REPLACE #include "../xml/
%       - FIND #include "dijkstra/      REPLACE #include "../dijkstra/
%   c) Build a Win32 or x64 Release of each project to compile the code.
%   d) Copy the generated .lib files to the following folder:
%
%   OPTI/Solvers/scip/Source/lib/win32 or win64
%
%   Next we need to copy all the required header files. Copy the header files
%   (and folders) of the following folders:
%       - blockmemshell
%       - nlpi
%       - objscip
%       - scip
%       - cppad (just configure.hpp)
%   to 
%
%   OPTI/Solvers/scip/Source/Include

% 3) MEX Interface
% The MEX interface supplied with SCIP (at the time of this development)
% was quite basic, thus has been updated for use with OPTI. Therefore I
% suggest you use the version supplied with OPTI.

% 4) ASL Interface
% If you wish to be able to solve AMPL .nl models using this interface, you
% can enable this functionality below with "haveASL". However you will need
% to complete building the ASL library, as detailed in Utilities/File
% IO/opti_AMPL_Install.m. Also ensure reader_nl.c and .h are placed in
% scip/Source/AMPL/. Note a bug exists in v3.0.1, so ensure the latest
% release is used.

% 5) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the SCIP MEX file. Once you have completed all the
% above steps, simply run this file to compile SCIP! You MUST BE in the 
% base directory of OPTI!

clear scip

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();
% Dependency Paths
ipoptdir = '..\..\ipopt\Source\';
mumpsdir = '..\..\mumps\Source\';
asldir = '..\..\..\Utilities\File IO\Source\';

% Switch to Enable AMPL ASL linking
haveASL = true;

fprintf('\n------------------------------------------------\n');
fprintf('SCIP MEX FILE INSTALL\n\n');
   
%Get IPOPT Libraries
post = [' -IInclude -I..\..\ipopt\Source\Include\Common\ -L' ipoptdir libdir ' -llibIpopt -llibmwma57']; 
%Get ASL libraries
if(haveASL)
   post = [post ' -I"' asldir 'Include" -L"' asldir libdir(1:end-1) '" -llibasl -DHAVE_ASL -DNO_STDIO1'];
end
%Get MUMPS Libraries
post = [post ' -I' mumpsdir 'Include -L' mumpsdir libdir '-llibdmumps_c -llibdmumps_f -llibseq_c -llibseq_f -llibmetis -llibpord'];
%Get Intel Fortran Libraries (for MUMPS build) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link];
%Get SCIP Includes and Libraries
post = [post ' -IInclude -IInclude/nlpi -IInclude/blockmemshell -L' libdir ' -llibscip -llibsoplex -llibut -output scip'];
   
%CD to Source Directory
cdir = cd;
cd 'Solvers/scip/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims scipmex.cpp scipeventmex.cpp scipnlmex.cpp';
if(haveASL)
    pre = [pre ' AMPL/reader_nl.c']; %include Stefan's ASL reader   
end
try
    eval([pre post])
    movefile(['scip.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:scip','Error Compiling SCIP!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');
