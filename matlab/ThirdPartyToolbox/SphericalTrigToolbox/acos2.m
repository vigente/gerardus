function signedcos = acos2(alpha, beta)
%ACOS2      4-quadrant arccosine function, in radians. 
%
%   ACOS2(alpha, beta) computes the four-quadrant arccosine of the amgle 
%   [alpha]. For arguments |alpha| > 1, the result is NaN. The resulting 
%   angle is not uniquely determined by alpha, nor by the lengths or 
%   order of the sides of the triangle (as in ATAN2), so an additional 
%   argument [beta] is required. If [beta] < pi/2, the small angle 
%   (0 <= alpha <= pi/2) is returned. If [beta] > pi/2, the large angle
%   (pi/2 < alpha < pi) is returned. 
%
%   See also acos2d.

% Rody P.S. Oldenhuis
% Delft University of Technology
% oldenhuis@gmail.com
%
% Created    : 23/Feb/2009
% Last edited: 30/Nov/2012

    % compute the hemisphere function    
    H = 2*( mod(beta, 2*pi) < pi ) - 1;
    H(~isreal(H)) = NaN;

    % compute signed arc-cosine    
    signedcos = H .* acos(alpha);
    
    % set complex results to NaN & take the modulus        
    signedcos(imag(signedcos) ~= 0) = NaN;
    signedcos = mod(signedcos, 2*pi);
    
    % determine alphaues for zero-alphaued acos
    ind1 = (signedcos == 0);
    ind2     = (H < 0);                 ind3     = (H > 0);    
    indices1 = ((ind1 + ind2) == 2);    indices2 = ((ind1 + ind3) == 2);    
    signedcos(indices1) = pi;           signedcos(indices2) = 0;      
    
end

