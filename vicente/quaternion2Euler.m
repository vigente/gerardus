%Function to minimize in order to convert a quaternion to Euler vectors.
%Actually, I guess there is a simple analytic solution to perform this
%conversion, but while I look for it, I can use a mean square difference
%minimization

function err = quaternion2Euler(angles, quatini)

    A=angles(1); B=angles(2); C=angles(3);
    %First rotation:
    qA = quaternion_axis_angle([0 0 0 1], A);
    
    %Second rotation:
    X1 = quaternion_product(qA, [0 1 0 0], quaternion_conj(qA));  %First calculate X'
    qB = quaternion_axis_angle(X1, B);
    
    %Third rotation:
    Z1 = quaternion_product(qA, [0 0 0 1], quaternion_conj(qA));
    Z2 = quaternion_product(qB, Z1, quaternion_conj(qB));
    qC = quaternion_axis_angle(Z2, C);
    
    quat = quaternion_product(qC,qB,qA);
    
    err = quatini - quat';