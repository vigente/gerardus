%%*******************************************************************
%%  Convert problem from SDPT3 format to DSDP format
%%
%%  [AC,b] = readsdpt3(blk,A,C,b)
%%
%%  Input: blk, A, C, b = data in SDPT3 format.
%%
%% 
%% DSDP: version 5.0
%% Copyright (c) 2004 by
%% S. Benson Y. Ye
%% Last modified: 2 Jan 04
%%******************************************************************

function [AC,b]=readsdpt3(blk,A,C,b)

AC=blk;
[p,m]=size(blk);
for bb=1:p
  conetype=blk{bb,1};
  if conetype=='s',
    [pp,mm]=size(A{bb});
    scl=ones(pp,1)/sqrt(2.0);
    dim=blk{bb,2};
    k=0;na=0;
    cc=sparse([]);
    for blki=1:length(dim)
      n=dim(blki);
      nb=na+n;
      for j=1:n,  k=k+j; scl(k)=1.0; end;
      cc=[cc;dvec(C{bb}(na+1:nb,na+1:nb))];
      na=nb;
    end;
    scl=spdiags(scl,0,k,k);
    AC{bb,3}=[scl*A{bb} cc];
    AC{bb,1}='SDP';
  elseif conetype=='l',
    AC{bb,3}=[A{bb} C{bb} ];
    AC{bb,1}='LP';
  elseif conetype=='q',
    error("Cannot do SOCP Constraints");
  end
end;
