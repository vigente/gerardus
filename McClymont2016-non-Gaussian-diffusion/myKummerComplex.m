function [chg]=myKummerComplex(a,b,z)

%
%KUMMERCOMPLEX   Confluent hypergeometric function 1F1 
% (Kummer function).
%
% KUMMERCOMPLEX(a,b,z) is the confluent hypergeometric 
% function 1F1 (Kummer function) for complex 
% parameters a, b and complex variable z.
%
% Example: 
%    KUMMERCOMPLEX(1+0.5i,2-3.1i,4+2i) 
%     equals   0.33300865268261 - 0.02369621687656i



%
% In general case the program calculates the sum of 
% convergent series defining the function until the next 
% term becomes too small (in comparison with the sum of all
% previous terms). The case of large abs(z) is considered 
% separately (e.g., see 13.5.1 in Abramowitz, Stegun 
% "Handbook of special functions", 1964). Some simple cases
% of integer parameter values are considered separately as 
% well.
%
% The function controls the loss of precision and makes 
% check for insufficient number of members in the series. 
% It prints warning if there are any problems. Otherwise, 
% if everything is ok,  the results seem to coincide with 
% Matematica 4.1 with 10-digit precision.
%
% This function is largely based at "Fortran library of 
% special functions" which was converted to Matlab. 
% Unfortunatey, the library can compute confluent 
% hypergeometric function only for real values of a and b. 
% So this file may be considered as its generalization 
% for complex a and b.
%
% This function also requires cgama.m file which computes 
% Gamma function for complex variables. This file was taken
% from just the same "Fortran library" and insignificantly 
% modified.
%


%% Darryl 5 Feb 2016 - this is all this function does for most voxels. Hopefully it is generalisable.
j = 1:5000;
crg = bsxfun(@times, (a+j-1)./(j.*(b+j-1)), z(:));
crg_prod = cumprod(crg,2);
chg = 1 + sum(crg_prod,2);
            

