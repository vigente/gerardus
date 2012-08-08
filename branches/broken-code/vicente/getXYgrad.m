function XYgrad = getXYgrad(phi)

XYgrad=zeros(size(phi,1), size(phi,2));

for i=2:size(phi,1)-1
		for j=2:size(phi,2)-1
			XYgrad(i,j)=(phi(i-1,j+1)+phi(i+1,j-1)-phi(i-1,j-1)-phi(i+1,j+1))/4;
        end
    end