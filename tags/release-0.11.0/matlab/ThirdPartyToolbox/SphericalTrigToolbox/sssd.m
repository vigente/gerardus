
function [A1, B1, C1, A2, B2, C2] = sssd(a, b, c)
%SSS   gives both solutions to the side-side-side problem, in degrees.
%
%   SSSD(a, b, c) may result in a vector filled with NaNs if the existence
%   condition |180 - a| - |180 - b| <= |180 - c| <= |180 - a| + |180 - b|
%   is not met. 
%
%   See also SSS.

% Rody P.S. Oldenhuis
% Delft University of Technology
% oldenhuis@gmail.com
%
% Crated     : 23/Feb/2009
% Last edited: 30/Nov/2012

    % find solutions by calling sss
    r2d = 180/pi;   
    d2r = 1/r2d;    
    [A1, B1, C1, A2, B2, C2] = sss(a*d2r, b*d2r, c*d2r);
    [A1, B1, C1, A2, B2, C2] = deal(A1*r2d, B1*r2d, C1*r2d, A2*r2d, B2*r2d, C2*r2d);
    
end