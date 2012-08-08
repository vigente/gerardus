function quat = versor2quaternion(versor)

quat=zeros(1,4);

for comp=2:4
    quat(comp) = versor(comp-1);
end

sinangle2 = norm(versor);
cosangle2 = sqrt( 1 - sinangle2^2 );
quat(1) = cosangle2;