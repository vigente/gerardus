function YYgrad = getYYgrad(phi)

YYgrad=zeros(size(phi,1), size(phi,2));

for i=2:size(phi,1)-1
		for j=1:size(phi,2)
				YYgrad(i,j)=(2*phi(i,j)-phi(i-1,j)-phi(i+1,j))/4;
        end
    end
