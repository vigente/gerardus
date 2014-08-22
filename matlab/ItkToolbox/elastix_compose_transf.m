function t3 = elastix_compose_transf(t1, t2)
% elastix_compose_transf  Composition of two transforms in elastix format.
%
% elastix_compose_transf composes two transforms produced by elastix.
%
% TC = elastix_compose_transf(T1, T2)
%
%   T1, T2 are two transforms in the format produced by elastix (see help
%   elastix for details).
%
%   TC is the composed transform.
%
%   transformix(TC, IM) produces the same result in one step than
%   transformix(T2, transformix(T1, IM)).
%
%   Implemented transforms:
%
%     'SimilarityTransform': TC has the same rotation centre as T1.
%
% See also: elastix, transformix, elastix_transf_imcoord2.

% Author: Ramon Casero <rcasero@gmail.com>
% Copyright © 2014 University of Oxford
% Version: 0.1.0
% $Rev$
% $Date$
% 
% University of Oxford means the Chancellor, Masters and Scholars of
% the University of Oxford, having an administrative office at
% Wellington Square, Oxford OX1 2JD, UK. 
%
% This file is part of Gerardus.
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. The offer of this
% program under the terms of the License is subject to the License
% being interpreted in accordance with English Law and subject to any
% action against the University of Oxford being under the jurisdiction
% of the English Courts.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

% check arguments
narginchk(2, 2);
nargoutchk(0, 1);

if (~strcmp(t1.HowToCombineTransforms, 'Compose') ...
        || ~strcmp(t2.HowToCombineTransforms, 'Compose'))
    error('HowToCombineTransforms must be ''Compose''')
end

% combine transforms
switch (t1.Transform)
    
    case 'SimilarityTransform' % similarity transform
        
        if (~strcmp(t2.Transform, 'SimilarityTransform'))
            error('If t1 is a SimilarityTransform, then t2 must be a SimilarityTransform')
        end
        
        if ((length(t1.TransformParameters) ~= 4) ...
                || (length(t2.TransformParameters) ~= 4))
            error('SimilarityTransform must have 4 parameters')
        end
        
        % center transforms
        c1 = t1.CenterOfRotationPoint';
        t1 = center_similarity_transf(t1);
        t2 = center_similarity_transf(t2);
        
        % nomenclature
        s1 = t1.TransformParameters(1);
        s2 = t2.TransformParameters(1);
        theta1 = t1.TransformParameters(2);
        theta2 = t2.TransformParameters(2);
        R1 = [cos(theta1) -sin(theta1); sin(theta1) cos(theta1)];
        R2 = [cos(theta2) -sin(theta2); sin(theta2) cos(theta2)];
        d1 = t1.TransformParameters(3:4)';
        d2 = t2.TransformParameters(3:4)';
        
        % combine uncentered transforms
        sc = s1 * s2;
        thetac = theta1 + theta2;
        dc = s1 * R1 * d2 + d1;
        
        % recenter combined transform, so that it has the same center as
        % the first transform
        dc = dc + (sc * (R1 * R2) - eye(2)) * c1;
        
        % format combined output in elastix struct
        t3 = t1;
        t3.TransformParameters = [sc thetac dc'];
        t3.CenterOfRotationPoint = c1';
        
    otherwise
        
        error('Transform not implemented')
        
end

end

function t = center_similarity_transf(t)

% transformation nomenclature
s =     t.TransformParameters(1);
theta = t.TransformParameters(2);
d =     t.TransformParameters(3:4)';
c =     t.CenterOfRotationPoint';
R = [cos(theta) -sin(theta);...
    sin(theta) cos(theta)];

% center transform
t.TransformParameters(3:4) = (d + (eye(2) - s * R) * c)';
t.CenterOfRotationPoint(:) = 0;

end
