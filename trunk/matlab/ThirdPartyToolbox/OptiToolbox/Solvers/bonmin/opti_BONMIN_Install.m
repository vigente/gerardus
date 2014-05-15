%% BONMIN Install for OPTI Toolbox
% Copyright (C) 2012 Jonathan Currie (I2C2)

% NOTE - From OPTI v1.73 BONMIN is now dynamically linked against CPLEX
% v12.5.0. This step is optional, and requires the user to have CPLEX
% installed and licensed on their system, as well as the same version. Two
% mex files will be created, one with Cplex and one without, as Cplex must
% be present to load the Cplex version, even if not used.

% Switch to Enable Linking against CPLEX dll
haveCPLEX = true;
% Build Local Solver Libs (i.e. CLP, CBC, IPOPT supplied with BONMIN dist) 
% (Sometimes there may be dependency problems with newer versions)
buildLocals = false;

% This file will help you compile Basic Open-source Nonlinear Mixed INteger
% programming (BONMIN) for use with MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012
% - Intel Compiler XE (FORTRAN)
% - Intel Math Kernel Library

% To recompile you will need to get / do the following:

% 0) Complete Compilation as per OPTI instructions for CLP, CBC, MUMPS, and
% IPOPT, in that order.

% 1) Get BONMIN
% BONMIN is available from http://www.coin-or.org/Bonmin/. Download 
% the source. 

% 2) Compile BONMIN
% The easiest way to compile CLP is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

% %% Visual Studio Builder Commands
% % NOTE the versions of CBC, CGL, IPOPT etc supplied with BONMIN must be the
% % same as those we are linking against! Otherwise paste new source over the
% % top of the BONMMIN supplied source.
% path = 'full path to BONMIN here'; %e.g. 'C:\Solvers\BONMIN'
% n = 1;
% %BONMIN
% sdir = [path '\src'];
% hdrs = {[path '\..\Cgl\src'],[path '\..\Cbc\src'],[path '\..\Clp\src'],[path '\..\CoinUtils\src'],[path '\..\Ipopt\src'],[path '\..\Osi\src'], [path '\..\BuildTools\headers']};
% name = 'libbonmin';
% opts = [];
% opts.exPP = {'BONMIN_BUILD','IPOPT_BUILD','FUNNY_MA57_FINT','WIN64','COINHSL_HAS_MA57','_CRT_SECURE_NO_WARNINGS'};
% opts.exclude = {'BonCurvatureEstimator.cpp','BonCurvBranchingSolver.cpp'};
% opts.exFolder = {'Apps','Interfaces\Ampl','Interfaces\Filter','Algorithms\Ampl'};
% VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% if(haveCPLEX)    
%     [~,cplx_inc] = opti_FindCplex();
%     %Bonmin with CPLEX
%     VSPRJ(n).sdir = VSPRJ(n-1).sdir; VSPRJ(n).hdrs = VSPRJ(n-1).hdrs; VSPRJ(n).opts=VSPRJ(n-1).opts; 
%     VSPRJ(n).name = 'libbonmincplex';
%     VSPRJ(n).hdrs = [VSPRJ(n).hdrs cplx_inc];
%     VSPRJ(n).opts.exPP = [VSPRJ(n).opts.exPP 'COIN_HAS_CPX'];
%     n = n + 1;    
%     % OsiCpx (must add Osi folder manually to include)
%     sdir = [path '\..\Osi\src\OsiCpx'];
%     hdrs = {[path '\..\CoinUtils\src'],[path '\..\BuildTools\headers'],cplx_inc};
%     name = 'libosicpx';
%     opts = [];
%     opts.exPP = {'COINUTILS_BUILD','_CRT_SECURE_NO_WARNINGS'};          
%     VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
% end
% if(buildLocals)
%     clppath = 'full path to CLP here'; %e.g. 'C:\Solvers\CLP'
%     cbcpath = 'full path to CBC here'; %e.g. 'C:\Solvers\CBC'
%     ippath = 'full path to IPOPT here'; %e.g. 'C:\Solvers\IPOPT'
%     mumpspath = 'full path to MUMPS here'; %e.g. 'C:\Solvers\MUMPS'
%     %CLP
%     sdir = [clppath '\src'];
%     hdrs = {[clppath '\..\CoinUtils\src'], [clppath '\..\Osi\src'], [clppath '\..\BuildTools\headers']};
%     name = 'libclp';
%     opts = [];
%     opts.exPP = {'CLP_BUILD','_CRT_SECURE_NO_WARNINGS'};
%     opts.exclude = {'MyEventHandler.cpp','MyMessageHandler.cpp','unitTest.cpp','CbcOrClpParam.cpp',...
%                     'ClpCholeskyMumps.cpp','ClpCholeskyUfl.cpp','ClpCholeskyWssmp.cpp','ClpCholeskyWssmpKKT.cpp'};
%     opts.exFilter = {'Abc*','CoinAbc*'};
%     VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
%     % CoinUtils
%     sdir = [clppath '\..\CoinUtils\src'];
%     hdrs = {[clppath '\..\BuildTools\headers']};
%     name = 'libcoinutils';
%     opts = [];
%     opts.exPP = {'COINUTILS_BUILD','_CRT_SECURE_NO_WARNINGS'};          
%     VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
%     % Osi
%     sdir = [clppath '\..\Osi\src'];
%     hdrs = {[clppath '\..\CoinUtils\src'] [clppath '\..\BuildTools\headers']};
%     name = 'libosi';
%     opts = [];
%     opts.exPP = {'OSI_BUILD','_CRT_SECURE_NO_WARNINGS'};   
%     opts.exFolder = {'OsiCpx','OsiGlpk','OsiGrb','OsiMsk','OsiSpx','OsiXpr','OsiCommonTest'};
%     VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
%     % CBC
%     sdir = [cbcpath '\src'];
%     hdrs = {[cbcpath '\..\CoinUtils\src'], [cbcpath '\..\Cgl\src'], [cbcpath '\..\Clp\src'], [cbcpath '\..\Osi\src'], [cbcpath '\..\BuildTools\headers']};
%     name = 'libcbc';
%     opts = [];
%     opts.exPP = {'CBC_BUILD','COIN_FAST_CODE','CLP_FAST_CODE','USE_CBCCONFIG','COIN_NO_TEST_DUPLICATE','_CRT_SECURE_NO_WARNINGS'};
%     opts.exclude = {'unitTest.cpp','unitTestClp.cpp','CoinSolve.cpp','CbcGeneric.cpp','CbcBranchBase.cpp'};
%     VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
%     % CGL
%     sdir = [cbcpath '\..\Cgl\src'];
%     hdrs = {[cbcpath '\..\CoinUtils\src'], [cbcpath '\..\Clp\src\'], [cbcpath '\..\Osi\src'], [cbcpath '\..\BuildTools\headers']};
%     name = 'libcgl';
%     opts = [];
%     opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};          
%     VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;
%     % IPOPT (remember to copy updated IpMa57TSolverInterface files and define WIN64 for x64 build)
%     sdir = [ippath '\src'];
%     hdrs = {[mumpspath '\include'], [mumpspath '\libseq'], [ippath '\..\BuildTools\headers']};
%     name = 'libipopt';
%     opts = [];
%     opts.exPP = {'IPOPT_BUILD','FUNNY_MA57_FINT','WIN64','COINHSL_HAS_MA57','_CRT_SECURE_NO_WARNINGS'};
%     opts.exFolder = {'Apps','Algorithm\Inexact','contrib\LinearSolverLoader'};
%     opts.exclude = 'AmplTNLP.cpp';
%     VSPRJ(n).sdir = sdir; VSPRJ(n).hdrs = hdrs; VSPRJ(n).name=name; VSPRJ(n).opts=opts; n = n + 1;    
% end
% %Write the Solution File
% VS_WriteSol(VSPRJ)
% %%

% Once complete, you will have a directory called BONMIN/libbonmin. Open the
% Visual Studio 2012 solution file, then complete the following steps:
%   a) In libbonmin open Interfaces/config_default.h and comment the line 
%   #define COIN_HAS_ASL 1
%   b) In libosicpx under project properties -> C/C++ -> General add the
%   Osi\src\Osi folder as an additional include directory.
%   a) Build a Win32 or x64 Release of each project to compile the code.
%   b) Copy the generated .lib files to the following folder:
%
%   OPTI/Solvers/bonmin/Source/lib/win32 or win64
%
%   You will also need to copy header files from the following directories
%   (keeping the same folder structure)
%   following folder:
%       - Algorithms
%       - Algorithms/OaGenerators
%       - CbcBonmin
%       - Interfaces
%   to:
%   OPTI/Solvers/bonmin/Source/Include/

% 3) BONMIN MEX Interface
% The BONMIN MEX Interface is based primarily on Peter Carbonetto's MEX
% Interface to IPOPT, with a few small changes I have made to enable use 
% with BONMIN. Notable changes include changes to the options + extra class
% methods in matlabprogram.

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the BONMIN MEX file. Once you have completed all 
% the above steps, simply run this file to compile BONMIN! You MUST BE in 
% the base directory of OPTI!

clear bonmin bonminCplex

% Modify below function if it cannot find Intel MKL on your system.
[mkl_link,mkl_for_link] = opti_FindMKL();
% Get Arch Dependent Library Path
libdir = opti_GetLibPath();
% Dependency Paths
clpdir = '..\..\clp\Source\';
cbcdir = '..\..\cbc\Source\';
mumpsdir = '..\..\mumps\Source\';
ipoptdir = '..\..\ipopt\Source\';
%Optionally reset dependency paths if using local library builds
if(buildLocals)
    clpdir = '';
    cbcdir = '';
    ipoptdir = '';
end

fprintf('\n------------------------------------------------\n');
fprintf('BONMIN MEX FILE INSTALL\n\n');

%Get CBC & CGL Libraries
post = [' -I' cbcdir 'Include\Cbc -I' clpdir 'Include\Osi -I' cbcdir 'Include\Cgl '];
post = [post ' -L' cbcdir libdir ' -llibCbc -llibCgl -llibut'];
%Get CLP & COINUTILS libraries
post = [post ' -I' clpdir 'Include -I' clpdir 'Include\Clp -I' clpdir 'Include\Coin'];
post = [post ' -L' clpdir libdir ' -llibClp -llibOsi -llibCoinUtils'];
%Get MUMPS Libraries
post = [post ' -I' mumpsdir '\Include -L' mumpsdir libdir ' -llibdmumps_c -llibdmumps_f -llibseq_c -llibseq_f -llibmetis -llibpord'];
%Get Intel Fortran Libraries (for MUMPS build) & MKL Libraries (for BLAS)
post = [post mkl_link mkl_for_link];
%Get IPOPT libraries
post = [post ' -I' ipoptdir 'Include\Common -I' ipoptdir 'Include\Interfaces -I' ipoptdir 'Include\LinAlg'];
post = [post ' -L' ipoptdir libdir ' -llibIpopt -llibmwma57 -DhaveMA57'];
%Get BONMIN Includes
post = [post ' -IInclude -IInclude/Interfaces -IInclude/CbcBonmin -IInclude/Algorithms -IInclude/Algorithms/OaGenerators']; 
%Output Common
post = [post ' -L' libdir ' -DBONMIN_BUILD -DIPOPT_BUILD -llibbonmin -output bonmin'];

%CD to Source Directory
cdir = cd;
cd 'Solvers/bonmin/Source';

%Compile & Move
pre = ['mex -v -largeArrayDims bonminmex.cpp iterate.cpp options.cpp bonminoptions.cpp matlabinfo.cpp '...
       'callbackfunctions.cpp matlabexception.cpp matlabfunctionhandle.cpp matlabprogram.cpp matlabjournal.cpp sparsematrix.cpp'];
try
    eval([pre post])
    movefile(['bonmin.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:bonmin','Error Compiling BONMIN!\n%s',ME.message);
end

%BUILD CPLEX BONMIN VERSION
if(haveCPLEX)
    % Modify below function it it cannot find IBM ILOG CPLEX on your system
    [CPLX_link] = opti_FindCplex();
    post = [post CPLX_link ' -llibosicpx -DHAVE_CPLEX' ];
    % Include BONMIN build with CPLEX
    post = regexprep(post,' -llibbonmin -output bonmin',' -llibbonminCplex -output bonminCplex');
    
    %Compile & Move
    pre = ['mex -v -largeArrayDims bonminmex.cpp iterate.cpp options.cpp bonminoptions.cpp matlabinfo.cpp '...
           'callbackfunctions.cpp matlabexception.cpp matlabfunctionhandle.cpp matlabprogram.cpp matlabjournal.cpp sparsematrix.cpp'];
    try
        eval([pre post])
        movefile(['bonminCplex.' mexext],'../','f')
        fprintf('Done!\n');
    catch ME
        cd(cdir);
        error('opti:bonmin','Error Compiling BONMIN CPLEX Version!\n%s',ME.message);
    end
end

cd(cdir);
fprintf('------------------------------------------------\n');
