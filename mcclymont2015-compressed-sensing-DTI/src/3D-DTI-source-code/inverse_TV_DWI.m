function [ res ] = inverse_TV_DWI( y )


sz = size(y);

if length(sz) == 4 % 2D
    res = zeros([sz(1), sz(2), sz(4)]);

    for i = 1:size(y, 4) % each dw volume
        res(:,:,i) = inverse_TV_2D(y(:,:,:,i));
    end
    
    
else % 3D
    res = zeros([sz(1), sz(2), sz(3), sz(5)]);

    for i = 1:size(y, 5) % each dw volume
        res(:,:,:,i) = adjDx(y(:,:,:,1,i)) + adjDy(y(:,:,:,2,i)) + adjDz(y(:,:,:,3,i));
    end
end
    
end


function res = adjDy(x)
    res = x(:,[1,1:end-1],:) - x;
    res(:,1,:) = -x(:,1,:);
    res(:,end,:) = x(:,end-1,:);
end

function res = adjDx(x)
    res = x([1,1:end-1],:,:) - x;
    res(1,:,:) = -x(1,:,:);
    res(end,:,:) = x(end-1,:,:);
end

function res = adjDz(x)
    res = x(:,:,[1,1:end-1]) - x;
    res(:,:,1) = -x(:,:,1);
    res(:,:,end) = x(:,:,end-1);
end
