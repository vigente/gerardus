function Ygrad = getYgrad(phi)

Ygrad=zeros(size(phi,1), size(phi,2));

	for i=2:size(phi,1)-1
		for j=1:size(phi,2)
				Ygrad(i,j)=(phi(i+1,j)-phi(i-1,j))/2;
        end
    end
