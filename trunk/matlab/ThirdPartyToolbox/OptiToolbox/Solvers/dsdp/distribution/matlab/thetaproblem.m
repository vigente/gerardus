%%*******************************************************
%% theta: Lovasz theta number. 
%%
%% (P)   min  Tr C*X
%%       s.t. X(i,j) = 0 if (i,j) is an edge of G, 
%%            Tr(X) = 1.                         
%%  b = e1, 
%%  C = -ones(n), 
%%  A1 = eye(n), Ak = ei*ej' + ej*ei', if (i,j) is an edge. 
%%-------------------------------------------------------
%%
%% [objval,X] = theta(G);
%%
%% G: adjacency matrix of a graph.
%%
%% DSDP: version 5.6 
%% Copyright (c) 2005 by
%% S.J. Benson and Y. Ye
%% Last modified: 3 Feb 05
%%*******************************************************

 function [objval,X] = thetaproblem(G);

    if ~isreal(G); error('only real G allowed');  end; 

    n = length(G); 
    nn=n*(n+1)/2;

    [idx,idy]=find(triu(G,1));
    nedges=length(idx);
    b = [1 zeros(1,nedges)]'; 
    AC=cell(1,3);
    AC{1,1} = 'SDP';  AC{1,2} = n; 
    AC{1,3} = sparse(nn,nedges+2);
    AC{1,3}(:,1) = dvec(speye(n)); 
    AC{1,3}(:,nedges+2) = -ones(nn,1);


    for k = 1:nedges 
       AC{1,3}(:,k+1) = dsparse([idx(k)],[idy(k)],[1],n,n);    
    end;

    y0 = -2*n*b;

    OPTIONS=doptions;
    OPTIONS.print=5;
    OPTIONS.r0 = 0;
    [STAT,y,Xv] = dsdp(b,AC,OPTIONS,y0);
    objval = dot(b,y);
    X=dmat(Xv{1});
%%=======================================================

