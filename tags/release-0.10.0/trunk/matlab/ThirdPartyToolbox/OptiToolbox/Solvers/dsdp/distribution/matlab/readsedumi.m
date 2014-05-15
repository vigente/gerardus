%%*******************************************************************
%%  Convert problem from SeDuMi format to DSDP format
%%
%%  [AC,b] = readsedumi(At,b,C,K)
%%
%%  Input: At, b, C, K = Data in SeDuMi format.
%%
%% 
%% DSDP: version 5.0
%% Copyright (c) 2004 by
%% S. Benson Y. Ye
%% Last modified: 2 Jan 04
%%******************************************************************

  function [AC,b]=readsedumi(At,bb,c,K);
%%
%%  First, load the matlab file containing At, c, b, and K
%%
  if (size(c,1) == 1), c = c'; end;
  if (size(bb,1) == 1), bb = bb'; end;
  [nn,mm] = size(At); if (max(size(c)) == 1); c = c*ones(nn,1); end; 

  if ~isfield(K,'l'); K.l = 0; end   
  if ~isfield(K,'q'); K.q = 0; end
  if ~isfield(K,'s'); K.s = 0; end
  if K.l == 0 | isempty(K.l); K.l = 0; end;
  if sum(K.q) == 0 | isempty(K.q); K.q = 0; end
  if sum(K.s) == 0 | isempty(K.s); K.s = 0; end 
%%
%%
%%
  AC=cell(1,3);
block=1;
top=1;

   m = length(bb);
   b = -bb;

   if (K.l > 0) 
      AC{block,1} = 'LP';
      AC{block,2} = K.l;
      A = At(top:top+n-1,:);
      CC = c(top:top+n-1,:);
      AC{block,3} = [-A(:,1:m) CC(:,1)];
      block = block+1;
      top = top+K.l;
   end
   if (K.q > 0) 
      error(' Cannot accept SOCP cones\n');
   end

   if (K.s > 0) 
      for i = 1:length(K.s)
         n = K.s(i);
         AC{block,1} = 'SDP';
         AC{block,2} = n;
         A = At(top:top+n^2-1,:);
         CC = c(top:top+n^2-1,:);
         indicies = triu(reshape(1:n^2,n,n));
         indicies = indicies(find(indicies));
         AC{block,3} = [-A(indicies,1:m) CC(indicies,1)];
         block = block+1;
         top = top+n*n;
       end
   end  
 
%%

