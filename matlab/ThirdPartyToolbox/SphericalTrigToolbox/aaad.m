function [a1, b1, c1, a2, b2, c2] = aaad(A, B, C)
%AAAD  gives both solutions to the angle-angle-angle problem, in degrees.
%
%   AAAD(A, B, C) will result in NaNs if the existence condition 
%   |pi - |A|-|B|| <= |C| <= pi - ||A| - |B||
%   is not met. 
%
%   See also AAAD.

% Rody P.S. Oldenhuis
% Delft University of Technology
% oldenhuis@gmail.com
%
% Crated     : 23/Feb/2009
% Last edited: 30/Nov/2012
    
    % find both solutions by calling aaa directly
    r2d = 180/pi;   
    d2r = 1/r2d;    
    [a1, b1, c1, a2, b2, c2] = aaa(A*d2r, B*d2r, C*d2r);
    [a1, b1, c1, a2, b2, c2] = deal(a1*r2d, b1*r2d, c1*r2d, a2*r2d, b2*r2d, c2*r2d);
 
end