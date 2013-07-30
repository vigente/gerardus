function [c1, A1, B1, c2, A2, B2] = sas(a, C, b)
%SAS   gives both solutions to the side-angle-side problem, in radians.
%
%   SAS(a, C, b) returns the remaining unknowns of the spherical triangle, 
%   [c1, A1, B1, c2, A2, B2]. 
%
%   See also SASD, ACOS2.  

% Rody P.S. Oldenhuis
% Delft University of Technology
% oldenhuis@gmail.com
%
% Created    : 23/Feb/2009
% Last edited: 30/Nov/2012

    % first solution
    c1 = acos2( cos(a).*cos(b) + sin(a).*sin(b).*cos(C),  C );
    A1 = acos2( (cos(a) - cos(b).*cos(c1))./(sin(b).*sin(c1)),  a );
    B1 = acos2( (cos(b) - cos(a).*cos(c1))./(sin(a).*sin(c1)),  b );
    
    % second solution
    c2 = 2*pi - c1; 
    A2 = mod(A1 + pi, 2*pi);
    B2 = mod(B1 + pi, 2*pi);

end    