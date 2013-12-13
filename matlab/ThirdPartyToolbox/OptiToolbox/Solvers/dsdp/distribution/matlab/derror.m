%%*************************************************************************
%%  Verify solution of DSDP solver by computing objective values and errors
%%
%%  [pobj,dobj,err] = derror(STAT,y,X,b,AC)
%%
%%  Input:  AC   = DSDP Constraint data
%%          b    = DSDP objective
%%          X    = DSDP solution matrices
%%          y    = DSDP solution vector
%%          STAT = DSDP solver status
%%
%%  Output: objective values and errors 
%%          err(1): || A * X - b ||_2 / (1+||b||_1)
%%          err(2): 0 if X is PSD, else -smallest eigevalue /(1+||b||_1)
%%          err(3): || C - A'y - S||_F /  (1+||C||_1)
%%          err(4): 0 if S is PSD, else -smallest eigevalue /(1+||c||_1)
%%          err(5): duality gap/(1+abs(primal obj)+abs(dual obj))
%%          err(6): tr(XS)/(1+abs(primal obj)+abs(dual obj))
%%
%% DSDP5.0 
%% Copyright (c) 2003 by
%% S. Benson and Y. Ye
%% Last modified: December 2003
%%******************************************************************

function [pobj,dobj,err]=derror(STAT,y,X,b,AC); 

[ncones,p]=size(AC);
m=length(b);
if (p~=3) 
   error('Dimension of first cell array not correct.'); 
end;

berr=[-b;0];
seigerr=0;
xeigerr=0;
tracexs=0;
for k=1:ncones,
   name=AC{k,1};
   if (length(name)==3 & name=='SDP'),
      pp=length(AC{k,2});
      tnz=0; 
      for jj=1:pp,
         j=AC{k,2}(jj);
         k1=tnz+1;  k2=tnz+j*(j+1)/2;
         AACC=AC{k,3}(k1:k2,:);
         XX=X{k}(k1:k2);
         berr=berr+dadotx(AACC,XX);
         svec=AACC*[-y;1];
         s=dmat(svec);
         seigmin=min(eig(s));
         seigerr=min(seigerr,seigmin);
         c=dmat(AACC(:,m+1));
         cnorm=norm(full(c));
         x=dmat(XX);
         xeigmin=min(eig(x));
         xeigerr=min(xeigerr,xeigmin);
         tracexs=tracexs+sum(sum(x.*s));
         tnz=tnz+j*(j+1)/2;
      end;
   elseif (length(name)==2 & name=='LP'),
     berr=berr+AC{k,3}'*X{k};
     svec=AC{k,3}*[-y; 1]; 
     seigmin=min(svec);
     seigerr=min(seigerr,seigmin);
     xeigmin=min(X{k});
     xeigerr=min(xeigerr,xeigmin);
     tracexs=tracexs+dot(X{k},svec);
   elseif (length(name)==2 & name=='LB'),
     l=[AC{k,3}];
     xl=[X{k}];
     berr=berr - [xl; dot(l,xl)];
     tracexs=tracexs+dot(y-l,xl);
     xeigerr=min(xeigerr,min(xl));
     seigerr=min(xeigerr,min(y-l));
   elseif (length(name)==2 & name=='UB'),
     u=[AC{k,3}];
     xu=[X{k}];
     berr=berr + [xu; dot(u,xu)];
     tracexs=tracexs+dot(u-y,xu);
     xeigerr=min(xeigerr,min(xu));
     seigerr=min(xeigerr,min(u-y));
   else
     error('Cone type not recognized: %s \n',name); 
  end;
end;
dobj=dot(b,y);
pobj=berr(m+1);
berr(m+1)=0;
derr1 = norm(berr)/(1+norm(b,'inf'));
tracexs=tracexs/(1+abs(pobj)+abs(dobj));
dualitygap=(pobj-dobj)/(1+abs(pobj)+abs(dobj));
err=[derr1,xeigerr,0,seigerr,dualitygap,tracexs];


if (strcmp(STAT.stype,'Unbounded'))
normS=seigerr;
err(3)=normS/dobj;
err(4)=err4/dobj;
err(5)=0;
err(6)=0;
end;

if (strcmp(STAT.stype,'Infeasible'))
err(1)=norm(berr)/(abs(pobj)+1);
err(2)=xeigmin/(abs(pobj)+1);
err(5)=0;
err(6)=0;
end;

