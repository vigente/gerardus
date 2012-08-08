function quat=matrix2quaternion(mat)
    T = 1 + trace(mat);
    if(T>.00000001)
        S = .5 / sqrt(T);
        quat(1) = .25 / S;
        quat(2) = ( mat(3,2) - mat(2,3) ) * S;
        quat(3) = ( mat(1,3) - mat(3,1) ) * S;
        quat(4) = ( mat(2,1) - mat(1,2) ) * S;
        
    elseif(max(diag(mat)) == mat(1,1))
        S  = sqrt( 1.0 + mat(1,1) - mat(2,2) - mat(3,3) ) * 2;
        quat(2) = 0.25 * S;
        quat(3) = (mat(1,2) + mat(2,1) ) / S;
        quat(4) = (mat(1,3) + mat(3,1) ) / S;
        quat(1) = (mat(2,3) - mat(3,2) ) / S;

    elseif(max(diag(mat)) == mat(2,2))
        S  = sqrt( 1.0 + mat(2,2) - mat(1,1) - mat(3,3) ) * 2;
        quat(2) = (mat(1,2) + mat(2,1) ) / S;
        quat(3) = 0.25 * S;
        quat(4) = (mat(2,3) + mat(3,2) ) / S;
        quat(1) = (mat(1,3) - mat(3,1) ) / S;

    else
        S  = sqrt( 1.0 + mat(3,3) - mat(1,1) - mat(2,2) ) * 2;
        quat(2) = (mat(1,3) + mat(3,1) ) / S;
        quat(3) = (mat(2,3) + mat(3,2) ) / S;
        quat(4) = 0.25 * S;
        quat(1) = (mat(1,2) - mat(2,1) ) / S;   
    end
    
    %Check here if the program is correct (I am not sure that there are no
    %errors in the equations)
    matr=quaternion2matrix(quat);
    if(matr ~= mat)
        fprintf('Original matrix=');
        mat
        fprintf('Matrix obtained from the quaternion=');
        matr
    end