function [stop, converged] = ls_stopping(tol)
% Yin Zhang, April, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County
% Modified J.Currie AUT May 2013

global probData

Hist = probData.Hist;

stop = 0; converged = 0; probData.message = [];

iter = size(Hist,2);
trerror = Hist(1,iter);
if (trerror < tol) 
   stop = 1; converged = 1; 
   probData.message = 'Converged'; return; 
end;

small = 1.e-5; 
if (iter <= 2), return; end;
blowup = Hist(1:5,iter) > max(tol,min(Hist(1:5,:),[],2))/small;
if iter > 5
   blowup = blowup | Hist(1,iter) > 1000*min(Hist(1,:)) | ...
   all(Hist(1,iter-4:iter) > 1.01*(Hist(1,iter-5:iter-1)));
end
if any(blowup)
   stop = 1; converged = -1;
   probData.message = detectinf(Hist,tol);
end;

nstall = 5; if (iter <= nstall), return; end;
latest = iter-nstall:iter;
h = Hist(1:5,latest); hmean = mean(h,2)';

stall = false(5,1);
for i = 1:5
  stall(i) = all(abs(h(i,:)-hmean(i)) < 1.e-3*hmean(i));
  if(i > 1), stall(i) = stall(i) & hmean(i) > small; end;
end
if any(stall)
   if trerror < small
      stop = 1; converged = 2;
      probData.message = 'Likely converged but error > tolerance'; 
   else
      stop = 1; converged = -1;
      probData.message = detectinf(Hist,tol);
   end
end;


function message = detectinf(Hist,tol)
% 

% Yin Zhang, April, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

msgncov = 'Not Converged:';
msgpinf = [msgncov ' Primal Infeasible'];
msgdinf = [msgncov ' Dual Infeasible'];

minrp = min(Hist(2,:)+Hist(4,:));
minrd = min(Hist(3,:)); 

if minrp < tol, message = msgdinf; return; end;
if minrd < tol, message = msgpinf; return; end;

if minrp < 10*tol && minrd > sqrt(tol)
   message = [msgdinf '?']; return;
end;
if minrd < 10*tol && minrp > sqrt(tol)
   message = [msgpinf '?']; return;
end;

iter = size(Hist,2);
if Hist(6,iter) < -1.e+10 && Hist(7,iter) < 1.e+6
   message = [msgdinf '??']; return;
end;
if Hist(7,iter) > 1.e+10 && Hist(6,iter) > -1.e+6
   message = [msgpinf '??']; return;
end;
message = [msgncov ' Both Primal and Dual Infeasible?'];
