%Function to minimize in order to convert a quaternion to Euler vectors.
%Actually, I guess there is a simple analytic solution to perform this
%conversion, but while I look for it, I can use a mean square difference
%minimization

function rot_mat = Euler2matrix(ang)
a=ang(1); b=ang(2); c=ang(3);

% rot_mat = [cos(b)*cos(c) -cos(a)*sin(c)+sin(a)*sin(b)*cos(c) sin(a)*sin(c)+cos(a)*sin(b)*cos(c); ...
%             cos(b)*sin(c) cos(a)*cos(c)+sin(a)*sin(b)*sin(c) -sin(a)*cos(c)+cos(a)*sin(b)*sin(c); ...
%             -sin(b)         sin(a)*cos(b)                       cos(a)*cos(b) ];

% rot_mat = [cos(b)*cos(c) -cos(a)*sin(c)+sin(a)*sin(b)*cos(c) sin(a)*sin(c)+cos(a)*sin(b)*cos(c); ...
%             cos(b)*sin(c) cos(a)*cos(c)+sin(a)*sin(b)*sin(c) -sin(a)*cos(c)+cos(a)*sin(b)*sin(c); ...
%             -sin(b)         sin(a)*cos(b)                       cos(a)*cos(b) ];

        
a=ang(1); b=ang(2); c=ang(3);

rotz = [cos(c) -sin(c) 0; sin(c) cos(c) 0; 0 0 1];
roty = [cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)];
rotx = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)];

rot_mat = rotz*rotx*roty;
    
    
