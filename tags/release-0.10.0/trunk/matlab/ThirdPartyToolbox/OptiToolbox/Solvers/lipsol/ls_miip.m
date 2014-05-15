function [x,y,z,s,w,info] = ls_miip(A,b,c,ub,opts)
% MIIP        - Main program for MIIP algorithm.

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County
% Modified J.Currie AUT May 2013

global probData

verb = opts.verb;

if(verb), fprintf('\n<<<<< This is the LIPSOL MIIP algorithm >>>>>\n'); end

[m,n] = size(A); probData.nt = n + probData.nub; 
bnrm = norm(b); cnrm = norm(c); unrm = [];
if (probData.Ubounds_exist), unrm = norm(ub); end

if(probData.Dense_cols_exist)
    Aabs = abs(A(:,probData.ispars)); 
else
    Aabs = abs(A); 
end
ls_symbfct(Aabs*Aabs',opts); 
clear Aabs;

[x,y,z,s,w] = initpoint(A,b,c,ub,bnrm,cnrm);
[Rxz,Rsw,dgap] = complementy(x,z,s,w);
[Rb,Rc,Ru,rb,rc,ru] = feasibility(A,b,c,ub,x,y,z,s,w);
trerror = errornobj(x,y,w,b,rb,bnrm,c,rc,cnrm,ub,ru,unrm,dgap);

iter = 0; converged = 0; mzeros = sparse(m,1); nzeros = sparse(n,1);
if(verb), fprintf('\n  Residuals:   Primal     Dual    U-bounds    Gap     TR_error\n'); end
if(verb), fprintf('  ---------------------------------------------------------\n'); end

while (iter <= opts.maxiter && toc(opts.t0) <= opts.maxtime)%%% Loop begins %%%
    if(opts.monitor), ls_monitor(iter,rb,rc,ru,dgap); end
    %Iteration Printing
    if(verb), fprintf('  Iter %4i:  ', iter); end
    if(verb), fprintf('%8.2e  %8.2e  %8.2e  %8.2e  %8.2e\n',rb,rc,ru,dgap,trerror); end
    if iter > 0 
        %Check For Convergence
        [stop, converged] = ls_stopping(opts.tol);
        if(stop), break; end
    end

    [vmid,xn1,sn1] = getmisc(x,z,s,w);
    [P,U] = getpu(A,vmid);

    %Compute Search Direction
    [dx,dy,dz,ds,dw] = direction(A,P,U,Rb,Rc,Ru,Rxz,Rsw,vmid,xn1,sn1,z,w,0,bnrm,1);
    %Ratio Test
    [ap,ad] = ratiotest(x,z,s,w,dx,dz,ds,dw);

    if (opts.tau0*ap < 1 || opts.tau0*ad < 1)
        mu = centering(x,z,s,w,dx,dz,ds,dw,ap,ad,dgap,trerror);
        Rxz = dx.*dz; Rsw = ds.*dw;
        [dx2,dy2,dz2,ds2,dw2] = direction(A,P,U,mzeros,nzeros,nzeros,Rxz,Rsw,vmid,xn1,sn1,z,w,mu,bnrm,2);
        dx = dx + dx2; dy = dy + dy2; dz = dz + dz2;
        ds = ds + ds2; dw = dw + dw2;
        [ap,ad] = ratiotest(x,z,s,w,dx,dz,ds,dw);
    end

    %Update Iterates and Check Feasibility
    [Rxz,Rsw,dgap,x,y,z,s,w,opts.phi0] = update(ap,ad,x,y,z,s,w,dx,dy,dz,ds,dw,trerror,opts);    
    [Rb,Rc,Ru,rb,rc,ru] = feasibility(A,b,c,ub,x,y,z,s,w);
    [trerror,rrb,rrc,rru,rdgap,objp,objd] = errornobj(x,y,w,b,rb,bnrm,c,rc,cnrm,ub,ru,unrm,dgap);

    probData.Hist = [probData.Hist [trerror rrb rrc rru rdgap objp objd]'];
    iter = iter + 1;
end %%% Loop ends %%%

if(iter >= opts.maxiter)
    converged = -3;
    probData.message = 'Exceeded Maximum Iterations';
elseif(toc(opts.t0) > opts.maxtime)
    converged = -3.5;
    probData.message = 'Exceeded Maximum Time';
end

info(1) = converged;
info(2) = iter;
info(3) = trerror;


function [x,y,z,s,w] = initpoint(A,b,c,ub,bnrm,cnrm)
% INITPOINT   - Specifying the starting point. 
% Usage: [x,y,z,s,w] = initpoint(A,b,c,ubounds,bnrm,cnrm)
% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global probData

idense = probData.idense;
ispars = probData.ispars;

[m,n] = size(A); y = zeros(m,1);
pmin = max(bnrm/100, 100);
dmin = cnrm*.425; dmin = max(dmin, pmin/40);
pmin = min(pmin, 1.e+4); dmin = min(dmin, 1.e+3);

e = ones(n,1); 
if probData.Dense_cols_exist
    ns = length(ispars); nd = length(idense); 
    Ds = sparse(1:ns,1:ns,e(ispars),ns,ns,ns);
    P = A(:,ispars)*Ds*A(:,ispars)';
    Dd = sparse(1:nd,1:nd,sqrt(e(idense)),nd,nd,nd);
    U = full(A(:,idense)*Dd);
    x = A'*densol(P,U,full(b),1); 
else
    P = A*sparse(1:n,1:n,e,n,n,n)*A';
    rho = min(100,bnrm);
    x = A'*ls_linsys(P,full(b-rho*A*e),1) + rho*e; 
end
pmin = max(pmin, -min(x)); x = max(pmin,x); 
z = full((c+dmin).*(c > 0) + dmin*(-dmin < c & c <= 0) - c.*(c <= -dmin));

s = []; w = [];
if (probData.Ubounds_exist)
   s = spones(ub).*max(pmin, ub-x);
   w = spones(ub).*( dmin*(c > 0) + ...
       (dmin - c).*(-dmin < c & c <= 0) - 2*c.*(c <= -dmin) );
end


function [Rxz,Rsw,dgap] = complementy(x,z,s,w)
% COMPLEMENTY - Evaluate the complementarity vectors and gap.
%	Usage:  [Rxz,Rsw,dgap] = complementy(x,z,s,w)
% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global probData

Rxz = x.*z; Rsw = [];
if (probData.Ubounds_exist), Rsw = s.*w; end
dgap = full(sum([Rxz; Rsw]));


function [Rb,Rc,Ru,rb,rc,ru] = feasibility(A,b,c,ub,x,y,z,s,w)
% FEASIBILITY - Evaluate feasibility residual vectors and their norms.
%	Usage:  [Rb,Rc,Ru,rb,rc,ru] = feasibility(A,b,c,ubounds)
% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global probData

Rb = A*x - b; Rc = A'*y + z - c; Ru = [];
if (probData.Ubounds_exist)
   Rc = Rc - w;
   Ru = spones(s).*x + s - ub;
end
rb = norm(Rb); rc = norm(Rc); ru = 0;
if(probData.Ubounds_exist), ru = norm(Ru); end
if(any(isnan([rb rc ru]))), error(' NaN occured.'); end


function [trerror,rrb,rrc,rru,rdgap,objp,objd] = errornobj(x,y,w,b,rb,bnrm,c,rc,cnrm,ubounds,ru,unrm,dgap)
% RELERROR    - Calculate the total relative error and objective values.

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global probData

rrb = rb/max(1,bnrm); rrc = rc/max(1,cnrm); rru = 0;
objp = full(c'*x); objd = full(b'*y);
if (probData.Ubounds_exist) 
   rru = ru/(1+unrm); 
   objd = objd - w'*ubounds;
end
rdgap = dgap/probData.nt;
%rdgap = abs(objp-objd)/(1 + min(abs(objp),abs(objd)));
trerror = max([rrb rrc rru rdgap]);


function [vmid,xn1,sn1] = getmisc(x,z,s,w)
% GETMISC     - Compute three quantities.
%             Usage: [vmid,xn1,sn1] = getmisc

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global probData

xn1 = reciprocal(x); sn1 = [];
if (~probData.Ubounds_exist)
   vmid = reciprocal(z.*xn1);
else
   sn1 = reciprocal(s);
   vmid = reciprocal(z.*xn1 + w.*sn1);
end;
cap = 1.e+11; vmid = full(min(cap,vmid));


function [P,U] = getpu(A,vmid)
% GETMISC     - Compute Matrix P and U
%             Usage: [P,U] = getpu(A,vmid)

% Yin Zhang, April, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global probData

idense = probData.idense; 
ispars = probData.ispars;
n = size(A,2);

if probData.Dense_cols_exist
   ns = length(ispars); nd = length(idense); 
   Ds = sparse(1:ns,1:ns,vmid(ispars),ns,ns,ns);
   P = A(:,ispars)*Ds*A(:,ispars)';
   Dd = sparse(1:nd,1:nd,sqrt(vmid(idense)),nd,nd,nd);
   U = full(A(:,idense)*Dd);
else
   P = A*sparse(1:n,1:n,vmid,n,n,n)*A'; U = [];
end;


function Y = reciprocal(X)
% RECIPROCAL  - Invert the nonzero entries of a matrix elementwise.
%             Y = RECIPROCAL(X) has the same sparsity pattern as X
%	      (except possibly for underflow).

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

if issparse(X)
   [m, n]  = size(X);
   [i,j,Y] = find(X);
   Y = sparse(i,j,1./Y,m,n);
else
   Y = 1./X;
end
ceiling = 1.e+16; Y = min(ceiling,Y);


function [dx,dy,dz,ds,dw] =  direction(A,P,U,Rb,Rc,Ru,Rxz,Rsw,vmid,xn1,sn1,z,w,mu,bnrm,flag)
% DIRECTION   - Computing search directions

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global probData

if (mu ~= 0), Rxz = Rxz - mu; end
Rc = Rc - Rxz.*xn1;
if (probData.Ubounds_exist)
   if (mu ~= 0), Rsw = Rsw - mu; end
   Rc = Rc + (Rsw - Ru.*w).*sn1;
end
rhs = -(Rb + A*(vmid.*Rc));

if probData.Dense_cols_exist
   dy = densol(P,U,full(rhs),flag);
else
    if(flag==1)
        dy=ls_linsys(P,full(rhs));
    else
        dy=ls_linsys([],full(rhs));
    end
end
dx = vmid.*(A'*dy + Rc);
dz = -(z.*dx + Rxz).*xn1;
ds = []; dw = [];
if (probData.Ubounds_exist)
   ds = -(dx.*spones(w) + Ru);
   dw = -(w.*ds + Rsw).*sn1;
end

resp = Rb + A*dx;
if norm(resp) > 1.e-1*bnrm
   dx = dx - A'*ls_solaat(full(resp));
end


function x = densol(P,U,b,flag)
% DENSOL      - Solve linear system with dense columns 
%               using either the Sherman-Morrison formula
%               or a preconditioned conjugate gradient method

% Yin Zhang, April, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County
% Revised Oct 2003, CAAM, Rice University

global probData

if(isempty(probData.Sherman_OK))
    probData.Sherman_OK = 1; 
end

bnrm = norm(b); x = zeros(size(b));
tol1  = min(1.e-2, 1.e-7*(1+bnrm));
tol2  = 10*tol1; tolcg = tol1;

iter = size(probData.Hist,2);
if iter > 0 
   trerror = probData.Hist(1,iter);
   tolcg = min(trerror, tolcg);
end

if probData.Sherman_OK
   x = sherman(P,U,b,flag); 
   resid = norm(b - P*x - U*(U'*x));
   if (resid > tol1), probData.Sherman_OK = 0; end
   if (resid < tol2), return; end
end
if ~probData.Sherman_OK 
   x = pcg(P,U,b,tolcg,flag,x);
end


function [ap,ad] = ratiotest(x,z,s,w,dx,dz,ds,dw)
% RATIOTEST   - Ratio test.
%          Usage: [ap,ad] = ratiotest(dx,dz,ds,dw)

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global probData

% ----- ratio test -----
ap = -1/min([dx./x; -0.1]);
ad = -1/min([dz./z; -0.1]);
if (probData.Ubounds_exist)
   as = -1/min([ds(s~=0)./nonzeros(s); -0.1]);
   aw = -1/min([dw(w~=0)./nonzeros(w); -0.1]);
   ap = min(ap, as); ad = min(ad, aw);
end


function mu = centering(x,z,s,w,dx,dz,ds,dw,ap,ad,dgap,trerror)
% CENTERING   - Computing centering parameter mu.
%             Usage: mu = centering(dx,dz,ds,dw,sn1,wn1,dgap)

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global probData

newdgap = (x + min(1,ap)*dx)'*(z + min(1,ad)*dz);
if (probData.Ubounds_exist) 
   newdgap = newdgap + (s + min(1,ap)*ds)'*(w + min(1,ad)*dw); 
end
sigmak = (newdgap/dgap)^2;

sigmin = 0; sigmax = .208; % Don't ask why.
p = ceil(log10(trerror)); 
if (p < -2 && dgap < 1.e+3), sigmax = 10^(p+1); end;
sigmak = max(sigmin, min(sigmax, sigmak));
%fprintf('   sigmak = %g\n',sigmak);

mu = sigmak*dgap/probData.nt;



function [Rxz,Rsw,dgap,x,y,z,s,w,phi0] = update(ap,ad,x,y,z,s,w,dx,dy,dz,ds,dw,trerror,opts)
% UPDATE      - Update the iterates.
%            UPDATE uses backtracking to make iterates stay in an
%            infinity-neighborhood of the central path.

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global probData

tau = opts.tau0;
if(~probData.backed), tau = .9 + 0.1*opts.tau0; end;
k = ceil(log10(trerror)); 
if (k <= -5), tau = max(tau,1-10^k); end;

ap = min(tau*ap,1); ad = min(tau*ad,1); 
xc = x; yc = y; zc = z; sc = s; wc = w;

step = [1 .9975 .95 .90 .75 .50];
for k = 1:length(step)
   x = xc + step(k)*ap*dx;
   y = yc + step(k)*ad*dy;
   z = zc + step(k)*ad*dz;
   s = sc + step(k)*ap*ds;
   w = wc + step(k)*ad*dw;
   [Rxz,Rsw,dgap] = complementy(x,z,s,w);
   phi = probData.nt*full(min([Rxz; Rsw(Rsw~=0)]))/dgap;
   if (max(ap,ad) == 1) || (phi >= opts.phi0), break; end;
end;
phi0 = min(opts.phi0, phi); 
if(k > 1), probData.backed = 1; end;


function x = sherman(P,U,b,flag)
% CHERMAN    Using the Sherman-Morrison formula to solve
%                 (P + U*U')x = b
%            where P is sparse and U is low-rank.

% Yin Zhang, April, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global probData

if (flag == 1)
   nd = size(U,2); 
   PinvU = zeros(size(U));
   for i = 1:nd
       PinvU(:,i) = ls_linsys(P,U(:,i)); 
   end;
   probData.Mdense = eye(nd) + U'*PinvU; 
else
   P = [];
end;

x = ls_linsys(P,b); 
tmp = U*(probData.Mdense\(U'*x));
x = x - ls_linsys(P,tmp);


function [x, iter] = pcg(P,U,b,tolcg,flag,x)
% pcg       Solve Ax = b using preconditioned conjugate gradient,
%           where A = P + U*U', P is sparse and U is dense.
%           P is used as the preconditioner.

n = length(b); bnrm = norm(b); iter = 0; 
maxiter = min(250,.5*n);
if (nargin < 6)
   x = zeros(n,1); r = b; rnrm = bnrm;
else
   r = b - P*x;
   if ~isempty(U), r = r - U*(U'*x); end;
   rnrm = norm(r); 
end;
while (rnrm > tolcg) && (iter <= maxiter)
    if (flag == 1), z = ls_linsys(P, r); end;
    if (flag ~= 1), z = ls_linsys([],r); end;
    iter = iter + 1;
    if iter == 1
       p = z; rtz = r'*z;
    else 
       oldrtz = rtz; rtz = r'*z;
       beta = rtz / oldrtz;
       p = z + beta*p;
    end;
    Ap = P*p; 
    if ~isempty(U), Ap = Ap + U*(U'*p); end;
    alpha = rtz / (p'*Ap);
    x = x + alpha*p;
    r = r - alpha*Ap;
    rnrm = norm(r);
end;
if iter == maxiter, error('Too many iterations in PCG solver'); end;