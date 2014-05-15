%% NOMAD Install for OPTI Toolbox

% This file will help you compile NOMAD for use with MATLAB. 

% My build platform:
% - Windows 8 x64
% - Visual Studio 2012

% To recompile you will need to get / do the following:

% 1) Get NOMAD
% NOMAD is available from http://www.gerad.ca/NOMAD/PHP_Forms/Download.php.
% Complete the download form then download the latest version. Once you
% have installed NOMAD, locate the /src/ directory.

% 2) Compile NOMAD
% The easiest way to compile NOMAD is to use the Visual Studio Project
% Builder included with OPTI. Use the following commands, substituting the
% required path on your computer:

% %% Visual Studio Builder Commands
% path = 'full path to NOMAD here'; %e.g. 'C:\Solvers\NOMAD'
% sdir = [path '\src'];
% name = 'libnomad';
% opts.exPP = {'_CRT_SECURE_NO_WARNINGS'};
% opts.exclude = 'nomad.cpp';
% VS_WriteProj(sdir,name,[],opts)
% %%

% Once complete, you will have a directory called NOMAD\libnomad. Open the
% Visual Studio 2012 project file, then complete the following steps:
%   a) In defines.hpp change the following: 
%       - line 194 change MAX_DIMENSION to 5000.
%   b) Build a Win32 or x64 Release to compile the code.
%   c) Copy the generated .lib file to the following folder:
%
%   OPTI/Solvers/nomad/Source/lib/win32 or win64
%
%   And rename it to libnomad.lib. 
%
%   You will also need to copy all header files from NOMAD/src to 
%   the following folder:
%
%   OPTI/Solvers/nomad/Source/Include

% 3) NOMAD MEX Interface
% % The NOMAD MEX Interface is a simple MEX interface I wrote to use NOMAD.
% Note the MEX interface contains two versions, a GERAD version (for the
% authors of NOMAD) and an OPTI version (for easier integration into this
% toolbox). By default the GERAD version is built, however to build the
% OPTI version simply add the preprocessor "OPTI_VERSION" to the compiler.

% 4) Compile the MEX File
% The code below will automatically include all required libraries and
% directories to build the NOMAD MEX file. Once you have completed all the
% above steps, simply run this file to compile NOMAD! You MUST BE in the 
% base directory of OPTI!

clear nomad

% Get Arch Dependent Library Path
libdir = opti_GetLibPath();

fprintf('\n------------------------------------------------\n');
fprintf('NOMAD MEX FILE INSTALL\n\n');

%Get NOMAD Libraries
post = [' -IInclude -L' libdir ' -llibnomad -llibut -output nomad -DOPTI_VERSION'];

%CD to Source Directory
cdir = cd;
cd 'Solvers/nomad/Source';

%Compile & Move
pre = 'mex -v -largeArrayDims nomadmex.cpp';
try
    eval([pre post])
    movefile(['nomad.' mexext],'../','f')
    fprintf('Done!\n');
catch ME
    cd(cdir);
    error('opti:nomad','Error Compiling NOMAD!\n%s',ME.message);
end
cd(cdir);
fprintf('------------------------------------------------\n');
