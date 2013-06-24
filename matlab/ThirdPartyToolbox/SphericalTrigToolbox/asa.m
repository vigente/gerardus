function [C1, a1, b1, C2, a2, b2] = asa(A, B, c)
%ASA   gives both solutions to the angle-side-angle problem, in radians.
%
%   ASA(A, B, c) returns the missing values C, a, b. It uses the
%   four-quadrant arccosine function ACOS2 to determine these values.
%
%   See also ACOS2, ASAD.

% Rody P.S. Oldenhuis
% Delft University of Technology
% oldenhuis@gmail.com
%
% Created    : 23/Feb/2009
% Last edited: 30/Nov/2012
    
    % first solution 
    % NOTE: normal acos (in stead of acos2) is indeed correct.
    C1 = acos( -cos(A) .*cos(B) + sin(A).*sin(B).*cos(c));
    a1 = acos( (cos(A) + cos(B).*cos(C1)) ./ (sin(B).*sin(C1)));
    b1 = acos( (cos(B) + cos(A).*cos(C1)) ./ (sin(A).*sin(C1)));
    
    C1(imag(C1) ~= 0) = NaN;
    a1(imag(a1) ~= 0) = NaN;
    b1(imag(b1) ~= 0) = NaN;
    
    % second solution
    C2 = 2*pi - C1;
    a2 = mod(a1 + pi, 2*pi);
    b2 = mod(b1 + pi, 2*pi);

end
