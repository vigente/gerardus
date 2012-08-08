function ang = matrix2Euler(rot_mat)

a=asin(rot_mat(3,2));

K=cos(a);

if(abs(a) > .00005)
    b = atan2( -rot_mat(3,1)/K, rot_mat(3,3)/K );
    c = atan2( -rot_mat(1,2)/K, rot_mat(2,2)/K );
else
    c = 0;
    b = atan2( rot_mat(1,1), rot_mat(2,1) );
end

ang = [a b c];
    
    
