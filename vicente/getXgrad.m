function Xgrad = getXgrad(phi)

Xgrad=zeros(size(phi,1), size(phi,2));

	for i=1:size(phi,1)
		for j=2:size(phi,2)-1
				Xgrad(i,j)=(phi(i,j+1)-phi(i,j-1))/2;
        end
    end

