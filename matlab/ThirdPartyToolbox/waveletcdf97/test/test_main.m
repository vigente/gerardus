%% Test script for waveletcdf97

% Pascal Getreuer 2006

fprintf('\n*** Perfect Reconstruction Test ***\n');
fprintf(['In this test, a random signal X is transformed one stage, then\n',...
      'inverse transformed, and the result is compared to the original.\n',...
      'Analytically, the transform is exactly inverted--the scheme is said\n',...
      'to have "perfect reconstruction."  This test verifies that this\n',...
      'property holds to machine precision.\n\n']);
Runs = 100;                    % Number of runs to average
N = ceil(logspace(...          % Lengths to test
   log10(15),log10(300),15));
AvgError = zeros(size(N));

for k = 1:length(N)
   % Create random input matrices
   X = rand(N(k),1,Runs);
   % Forward transform followed by inverse
   R = waveletcdf97(waveletcdf97(X,1),-1);
   % Compute the average error
   AvgError(k) = mean(max(abs(permute(X - R,[1,3,2])),[],1));
   fprintf('  length(X) = %3d   Error = %.2e\n',N(k),AvgError(k));
end

fprintf('\n');

figure;
plot(N,AvgError,'.-');
xlabel('length(X)');
ylabel('Error (average of 100 runs)');
title('Perfect Reconstruction Test');

fprintf('\n*** Vanishing Moments Test ***\n');
fprintf(['The CDF 9/7 wavelet is designed such that where the input signal\n',...
      'is locally a polynomial of cubic degree or lower, the resulting\n',...
      'detail (highpass) coefficients are equal to zero.  A wavelet is said\n',...
      'to have "N vanishing moments" if it has this property on polynomials\n',...
      'up to degree N-1, so CDF 9/7 has 4 vanishing moments.  This test\n',...
      'transforms a piecewise polynomial signal and displays the largest\n',...
      'detail coefficient magnitudes, verifying that the vanishing moments\n',...
      'hold to reasonable accuracy.\n\n']);

N = 64;
t = (0:N-1)/N;
X = [t.^0,t.^1,t.^2,t.^3];
Y = waveletcdf97(X,1);

Moment0 = norm(Y(2*N+2:2.5*N-2),inf);   % Largest detail coefficient from t^0
Moment1 = norm(Y(2.5*N+2:3*N-2),inf);   % Largest detail coefficient from t^1
Moment2 = norm(Y(3*N+2:3.5*N-2),inf);   % Largest detail coefficient from t^2
Moment3 = norm(Y(3.5*N+2:4*N-2),inf);   % Largest detail coefficient from t^3

fprintf('  Locally constant    Largest detail coefficient = %.2e\n',Moment0);
fprintf('  Locally linear      Largest detail coefficient = %.2e\n',Moment1);
fprintf('  Locally quadratic   Largest detail coefficient = %.2e\n',Moment2);
fprintf('  Locally cubic       Largest detail coefficient = %.2e\n\n',Moment3);


figure;
subplot(2,1,1);
plot(X);
axis([1,4*N,0,1.5]);
title('Vanishing Moments Test');
ylabel('Input signal');
text(0.2*N,1.2,'Constant');
text(1.2*N,1.2,'Linear');
text(2.2*N,1.2,'Quadratic');
text(3.2*N,1.2,'Cubic');
set(gca,'XTick',[1,N,2*N,3*N,4*N],'YTick',[0,1]);


subplot(2,1,2);
plot(Y(2*N+1:4*N));
xlim([1,2*N]);
ylabel('Detail coefficients');
set(gca,'XTick',[1,0.5*N,N,1.5*N,2*N]);

ArrowX = [0,0,0.022,0,-0.022]*N;
ArrowY = [0.08,0.01,0.03,0.01,0.03];

text(0.1*N,0.12,sprintf('%.2e',Moment0));
text(0.6*N,0.12,sprintf('%.2e',Moment1));
text(1.1*N,0.12,sprintf('%.2e',Moment2));
text(1.6*N,0.12,sprintf('%.2e',Moment3));
set(line(0.25*N+ArrowX,ArrowY),'Color',[0,0,0]);
set(line(0.75*N+ArrowX,ArrowY),'Color',[0,0,0]);
set(line(1.25*N+ArrowX,ArrowY),'Color',[0,0,0]);
set(line(1.75*N+ArrowX,ArrowY),'Color',[0,0,0]);
