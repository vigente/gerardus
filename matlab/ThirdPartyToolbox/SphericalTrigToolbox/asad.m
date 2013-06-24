function [C1, a1, b1, C2, a2, b2] = asad(A, B, c)
%ASAD   gives both solutions to the angle-side-angle problem, in degrees.
%
%   ASAD(A, B, c) returns the missing values C, a, b. It uses the
%   four-quadrant arccosine function ACOS2D to determine these values.
%
%   See also ACOS2, ACOS2D, ASA

% Rody P.S. Oldenhuis
% Delft University of Technology
% oldenhuis@gmail.com
%
% Crated     : 23/Feb/2009
% Last edited: 30/Nov/2012

    % find both solutions by calling asa directly    
    r2d = 180/pi;   
    d2r = 1/r2d;    
    [C1, a1, b1, C2, a2, b2] = asa(A*d2r, B*d2r, c*d2r);
    [C1, a1, b1, C2, a2, b2] = deal(C1*r2d, a1*r2d, b1*r2d, C2*r2d, a2*r2d, b2*r2d);

end


