function [B1, C1, c1, B2, C2, c2] = ssad(a, b, A)
%SSAD  gives both solutions to the side-side-angle problem, in degrees.
%
%   SSAD(a, b, A) will result in NaNs if the existence condition 
%   |sin b * sin A| <= | sin a |  is not met. 
%
%   See also SSA, MALD, MSLD.

% Rody P.S. Oldenhuis
% Delft University of Technology
% oldenhuis@gmail.com
%
% Crated     : 23/Feb/2009
% Last edited: 30/Nov/2012
    
    % find both solutions by calling ssa directly
    r2d = 180/pi;   
    d2r = 1/r2d;    
    [B1, C1, c1, B2, C2, c2] = ssa(a*d2r, b*d2r, A*d2r);
    [B1, C1, c1, B2, C2, c2] = deal(B1*r2d, C1*r2d, c1*r2d, B2*r2d, C2*r2d, c2*r2d);
    
end