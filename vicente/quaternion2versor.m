function versor = quaternion2versor(quat)

if(abs(norm(quat) - 1) > .0001)
    fprintf('Error: non-unit quaternion');
else
    versor = quat(2:4);
end
