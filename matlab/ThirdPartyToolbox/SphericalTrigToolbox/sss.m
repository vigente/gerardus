function [A1, B1, C1, A2, B2, C2] = sss(a, b, c)
%SSS   gives both solutions to the side-side-side problem, in radians.
%
%   SSS(a, b, c) results in NaNs for those indices where the existence 
%   condition |pi - a| - |pi - b| <= |pi - c| <= |pi - a| + |pi -b| is not
%   met. 
%
%   See also SSSD.

% Rody P.S. Oldenhuis
% Delft University of Technology
% oldenhuis@gmail.com
%
% Created    : 23/Feb/2009
% Last edited: 30/Nov/2012

    % first solution
    A1 = acos2( (cos(a) - cos(b).*cos(c))./(sin(b).*sin(c)), a);
    B1 = acos2( (cos(b) - cos(a).*cos(c))./(sin(a).*sin(c)), b);
    C1 = acos2( (cos(c) - cos(a).*cos(b))./(sin(a).*sin(b)), c);
    
    % second solution
    A2 = 2*pi - A1;
    B2 = 2*pi - B1;
    C2 = 2*pi - C1;
    
    % check constraints
    indices = ( ...
        (abs(pi-a) - abs(pi-b)) > abs(pi-c) | ...
        abs(pi-c) > (abs(pi-a) + abs(pi-b)) );
    A1(indices) = NaN;   B1(indices) = NaN;   C1(indices) = NaN;
    A2(indices) = NaN;   B2(indices) = NaN;   C2(indices) = NaN;
    
end