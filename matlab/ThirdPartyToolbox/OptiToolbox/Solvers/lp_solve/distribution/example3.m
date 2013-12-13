X = 0:125;
Y1 = (15000 - 120.*X)./210;
area(X, Y1)
Y2 = max((4000 - 110.*X)./30, 0);
Y3 = max(75 - X, 0);
Ytop = min([Y1; Y2; Y3]);
area(X, Ytop)
area(X, Ytop); axis([0 40 0 75])
hold on
[U V] = meshgrid(0:40, 0:75);
contour(U, V, 143.*U + 60.*V); hold off
x = [1 1; 110 30] \ [75; 4000]
format bank
P = [143 60] * x