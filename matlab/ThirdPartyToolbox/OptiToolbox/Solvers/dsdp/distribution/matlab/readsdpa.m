%%*******************************************************************
%%  Read in a problem in SDPA sparse format.
%%
%%  [AC,b] = readsdpa(fname)
%%
%%  Input: fname = name of the file containing SDP data in
%%                 SDPA foramt. 
%%
%% DSDP5
%% Copyright (c) 2005 by
%% S.J. Benson, Y. Ye
%% Last modified: 4 Feb 05
%%******************************************************************

   function [AC,b] = readsdpa(fname); 

%%
%%  Open the file for input
%%
   compressed = 0; 
  dblocktol=0;

   if exist(fname)
      fid = fopen(fname,'r');
   elseif exist([fname,'.Z']); 
      compressed = 1; 
      unix(['uncompress ',fname,'.Z']);
      fid = fopen(fname,'r');
   elseif exist([fname,'.gz']); 
      compressed = 2; 
      unix(['gunzip ',fname,'.gz']);
      fid = fopen(fname,'r');
   else
      error('** Problem not found, please specify the correct path.');
   end;
%%
%%  Clean up special characters and comments from the file 
%%
   [tmpr,count] = fscanf(fid,'%c');
   linefeeds = findstr(tmpr,char(10));
   comment_chars = '*"=';
   cumidx = [];
   for i=1:length(comment_chars)
      idx = findstr(tmpr,comment_chars(i));
      cumidx = [cumidx,idx];
   end;
   for j=length(cumidx):-1:1
      if (cumidx(j)==1) | (strcmp(tmpr(cumidx(j)-1),char(10)))
         tmpr(cumidx(j):linefeeds(min(find(cumidx(j)<linefeeds))))='';
      else
         tmpr(cumidx(j):linefeeds(min(find(cumidx(j)<linefeeds)))-1)='';
      end;
   end;
   special_chars = ',{}()';
   cumidx=[];
   for i=1:length(special_chars)
      idx = findstr(tmpr,special_chars(i));
      cumidx = [cumidx,idx];
   end;
   tmpr(cumidx) = blanks(length(cumidx));
   nztotal = length(tmpr);
%%
%% Close the file 
%%
   fclose('all');
   if compressed==1; unix(['compress ',fname]); end;
   if compressed==2; unix(['gzip ',fname]); end; 
%% 
%%  Next, read in basic problem size parameters.
%%
   datavec = sscanf(tmpr,'%f'); clear tmpr;
   if size(datavec,1) < size(datavec,2); datavec = datavec'; end; 
   m = datavec(1); 
   numblk  = datavec(2);
   blksize = datavec(2+[1:numblk]); 
   if size(blksize,1) > size(blksize,2); blksize = blksize'; end
%%
%% Get input  b.
%%
   idxstrt = 2+numblk; 
   b = datavec(idxstrt+[1:m]);    
   idxstrt = idxstrt+m; 
   b = -b;
%%
%% Construct AC
%%
   deblkidx = find( (blksize<0) | (blksize > dblocktol) ); 
   denumblk = length(deblkidx); 

      AC=cell(denumblk,3);

   for p = 1:length(deblkidx)
      n = blksize(deblkidx(p)); 
      if (n > 0); 
        AC{p,1}='SDP'; AC{p,2}=n; AC{p,3}=sparse(n*(n+1)/2,m+1);
      else
        AC{p,1}='LP'; AC{p,2}=abs(n); AC{p,3}=sparse(abs(n),m+1);
      end  
   end
%%
%% Construct single blocks of A,C
%%
   len = length(datavec);    
   Y = reshape(datavec(idxstrt+1:len),5,(len-idxstrt)/5)';   
   clear datavec;    
   Y = sortrows(Y,[1 2]); 
   matidx = [0; find(diff(Y(:,1)) ~= 0); size(Y,1)];
%%
   for k = 1:length(matidx)-1
      idx = [matidx(k)+1 : matidx(k+1)];       
      Ytmp  = Y(idx,1:5); 
      matno = Ytmp(1,1); 
      if (matno == 0) matno=m+1; end;
      for p = 1:denumblk 
         n = blksize(deblkidx(p));   
         idx = find( Ytmp(:,2)==deblkidx(p) ); 
         Yp = Ytmp(idx,3:5); 
         len = length(idx); 
           if (n>0)
               AC{p,3}(:,matno)= -dsparse(Yp(:,1),Yp(:,2),Yp(:,3),n,n); 
            else
               tmp = -sparse(Yp(:,1),ones(len,1),Yp(:,3),abs(n),1);
               AC{p,3}(:,matno)=tmp;
            end
       end
   end 
%%

%% 
%%******************************************************************








