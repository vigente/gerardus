function XXgrad = getXXgrad(phi)

XXgrad=zeros(size(phi,1), size(phi,2));

	for i=1:size(phi,1)
		for j=2:size(phi,2)-1
				XXgrad(i,j)=(2*phi(i,j)-phi(i,j-1)-phi(i,j+1))/4;
        end
    end

