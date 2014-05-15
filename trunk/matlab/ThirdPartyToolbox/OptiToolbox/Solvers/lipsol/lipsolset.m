function options = lipsolset(varargin)
%LIPSOLSET  Create or alter the options for Optimization with LIPSOL
%
% options = lipsolset('param1',value1,'param2',value2,...) creates an LIPSOL
% options structure with the parameters 'param' set to their corresponding
% values in 'value'. Parameters not specified will be set to the LIPSOL
% default.
%
% options = lipsolset(oldopts,'param1',value1,...) creates a copy of the old
% options 'oldopts' and then fills in (or writes over) the parameters
% specified by 'param' and 'value'.
%
% options = lipsolset() creates an options structure with all fields set to
% lipsolset defaults.
%
% lipsolset() prints a list of all possible fields and their function.
%
% See supplied LIPSOL Documentation for further details of these options.

%   Copyright (C) 2013 Jonathan Currie (I2C2)

% Print out possible values of properties.
if((nargin == 0) && (nargout == 0))
    printfields();
    return
end

%Names and Defaults
Names = {'phi0';'tau0';'monitor';'verb';'big';'cache_size';'loop_level'};
Defaults = {1e-5,0.9995,0,1,1e32,16,8}';        

%Enter and check user args
try
    options = opticheckset(Names,Defaults,@checkfield,varargin{:});
catch ME
    throw(ME);
end


function checkfield(field,value)
%Check a field contains correct data type
switch lower(field)
    %Scalar double
    case 'big'
        err = opticheckval.checkScalarDbl(value,field);
    %Scalar non negative double
    case 'phi0'
        err = opticheckval.checkScalarNonNeg(value,field);  
    %Scalar double 0 < val <= 1
    case 'tau0'
        err = opticheckval.checkScalarBoundLLE(value,field,0,1);
    %Scalar 0/1
    case 'monitor'
        err = opticheckval.checkScalar01(value,field); 
    %Scalar Int 0-2
    case 'verb'
        err = opticheckval.checkScalarIntBoundLELE(value,field,0,2);
    %Scalar Int 1 - 100
    case 'cache_size'
        err = opticheckval.checkScalarIntBoundLELE(value,field,1,100);
    %Scalar Int 1-8
    case 'loop_level'
        err = opticheckval.checkScalarIntBoundLELE(value,field,1,8);
    otherwise  
        err = MException('OPTI:SetFieldError','Unrecognized parameter name ''%s''.', field);
end
if(~isempty(err)), throw(err); end


function printfields()
%Print out fields with defaults
fprintf('       phi0: [ Relative Termination Tolerance: {1e-5} ] \n');
fprintf('       tau0: [ Initial Line-Search Step-Size Multiplier: {0.9995} ] \n');
fprintf('        big: [ Bound Infinite Value: {1e32} ] \n');
fprintf('    monitor: [ Graphical Iteration Monitor: On (1), Off {0} ] \n');
fprintf('       verb: [ Verbosity: On (1), Off {0} ] \n');
fprintf(' cache_size: [ SYMFCT Cache Size: {16} ] \n');
fprintf(' loop_level: [ BLKFCT Loop Unroll Level: 1,2,4,{8} ] \n');

fprintf('\n');
