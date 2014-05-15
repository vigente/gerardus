function [b1, c1, C1, b2, c2, C2] = aas(A, B, a)
%AAS   gives both solutions to the angle-angle-side problem, in radians.
%
%   AAS(A, B, a) may result in a vector filled with NaNs if the existence
%   condition |sin(B)sin(a)| <= |sin(A)| is not met. This function uses the
%   Middle Side Law function MSL.m and Middle Angle Law function MAL.m to
%   determine the solutions. 
%
%   See also MSL, MAL, AASD.

% Rody P.S. Oldenhuis
% Delft University of Technology
% oldenhuis@gmail.com
%
% Created    : 23/Feb/2009
% Last edited: 30/Nov/2012

    % first solution
    b0 = asin( (sin(B).*sin(a))./sin(A) );
    b0(imag(b0) ~= 0) = NaN;
    
    b1 = mod(b0, 2*pi);
    c1 = msl(a, b1, A, B);
    C1 = mal(A, B, a, b1);

    % second solution
    b2 = mod(pi - b1, 2*pi);
    c2 = msl(a, b2, A, B);
    C2 = mal(A, B, a, b2);
    
    % check constraints
    indices = ( abs(sin(B).*sin(a)) > abs(sin(A)) );
    b1(indices) = NaN;    c1(indices) = NaN; 
    C1(indices) = NaN;    b2(indices) = NaN; 
    c2(indices) = NaN;    C2(indices) = NaN;    

end


% Middle-angle-law
function C = mal(A, B, a, b)
%MAL    Computes the missing angle in a spherical triangle, in radians.
%
%   MAL(A, B, a, b) is the implementation of the Middle Angle Law, and
%   returns the missing angle C. 
%
%   See also MSL, MALD.

% Rody P.S. Oldenhuis
% Delft University of Technology
% Last edited: 23/Feb/2009

    % sine & cosine of C
    % NOTE: denomenator not needed
    sinC  =  sin(A).*cos(B).*cos(b) + sin(B).*cos(A).*cos(a);    
    cosC  = -cos(A).*cos(B) + sin(A).*sin(B).*cos(a).*cos(b);

    % C is the arctangent of the ratio of these two
    C = mod( atan2(sinC, cosC), 2*pi);

end

% Middle-side-law
function c = msl(a, b, A, B)
%MSL    Computes the missing side in a spherical triangle, in radians.
%
%   MSL(a, b, A, B) is the implementation of the Middle Side Law, and
%   returns the missing angular side c. 
%
%   See also MAL, MSLD.

% Rody P.S. Oldenhuis
% Delft University of Technology
% Last edited: 23/Feb/2009
    
    % sine & cosine c
    % NOTE: divisor not needed
    sinc  = (sin(a).*cos(b).*cos(B) + sin(b).*cos(a).*cos(A));	    
    cosc  = (cos(a).*cos(b) - sin(a).*sin(b).*cos(A).*cos(B));                

    % c is the arctangent of the sine over the cosine
    c = mod( atan2(sinc, cosc), 2*pi);

end

    