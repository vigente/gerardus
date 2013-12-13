%%*******************************************************************
%%  Read in a problem in SeDuMi format and convert to DSDP format
%%
%%  [AC,b] = readsedumifile(fname)
%%
%%  Input: fname.mat = name of the file containing SDP data in
%%                     SeDuMi format.
%%
%%  Input: AC  = DSDP formatted block structure and constraint data
%%         b   = Dual objective vector
%%
%% DSDP: version 5.0
%% Copyright (c) 2004 by
%% S. Benson Y. Ye
%% Last modified: 12 Jan 04
%%******************************************************************

  function [AC,b]=readsedumifile(fname);
%%
%%  First, load the matlab file containing At, c, b, and K
%%
  K.q = [];
  K.l = [];
  A = 0;
  At = 0;
  compressed = 0; 
  if exist([fname,'.mat.gz']); 
     compressed = 1; 
     unix(['gunzip ', fname,'.mat.gz']);
  elseif exist([fname,'.gz']); 
     compressed = 2; 
     unix(['gunzip ', fname,'.gz']);
  elseif exist([fname,'.mat.Z']); 
     compressed = 3; 
     unix(['uncompress ', fname,'.mat.Z']);
  elseif exist([fname,'.Z']); 
     compressed = 4; 
     unix(['uncompress ', fname,'.Z']);
  end
  if exist([fname,'.mat']) | exist(fname) 
     eval(['load ', fname]);
  else
     error('** Problem not found, please specify the correct path.');
  end

  if (compressed == 1)
     unix(['gzip ', fname,'.mat']);
  elseif (compressed == 2)
     unix(['gzip ', fname]);
  elseif (compressed == 3)
     unix(['compress ', fname,'.mat']);
  elseif (compressed == 4)
     unix(['compress ', fname]);
  end
%%
  if (size(c,1) == 1), c = c'; end;
  if (size(b,1) == 1), b = b'; end;
  if (At == 0), At = A'; end;  clear A;
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
%%  Then convert formats

  [AC b]=readsedumi(At,b,c,K);

return;
