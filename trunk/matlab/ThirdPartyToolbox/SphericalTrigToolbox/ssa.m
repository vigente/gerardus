function [B1, C1, c1, B2, C2, c2] = ssa(a, b, A)
%SSA   gives both solutions to the side-side-angle problem, in radians.
%
%   SSA(a, b, A) will result in NaNs if the existence condition 
%   |sin b * sin A| <= | sin a |  is not met. 
%
%   See also SSAD, MAL, MSL.

% Rody P.S. Oldenhuis
% Delft University of Technology
% oldenhuis@gmail.com
%
% Created    : 23/Feb/2009
% Last edited: 30/Nov/2012

    % first solution
    B0 = asin(sin(b).*sin(A)./sin(a));
    B0(imag(B0) ~= 0) = NaN;
    
    B1 = mod(B0, 2*pi);
    C1 = mal(A, B1, a, b);
    c1 = msl(a, b, A, B1);
    
    % second solution
    B2 = mod(pi - B1, 2*pi);    
    C2 = mal(A, B2, a, b);
    c2 = msl(a, b, A, B2);
    
    % check constraints    
    indices = ( abs(sin(b).*sin(A)) > abs(sin(a)) );
    B1(indices) = NaN;    C1(indices) = NaN;
    c1(indices) = NaN;    B2(indices) = NaN;
    C2(indices) = NaN;    c2(indices) = NaN;
    
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
    sinC  = sin(A).*cos(B).*cos(b) + sin(B).*cos(A).*cos(a);    
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