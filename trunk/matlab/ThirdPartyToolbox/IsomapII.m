function [Y, R, E] = IsomapII(D, n_fcn, n_size, options); 
% ISOMAPII   Computes Isomap embedding using an advanced version of
%             the algorithm in Tenenbaum, de Silva, and Langford (2000), 
%             which can take advantage of sparsity in the graph and 
%             redundancy in the distances. 
%
% [Y, R, E] = isomapII(D, n_fcn, n_size, options); 
%
% Input:
%    D = input-space distances between pairs of N points, which can 
%     take 1 of 3 forms: 
%      (1) a full N x N matrix (as in isomap.m)  
%      (2) a sparse N x N matrix (missing entries are treated as INF)
%      (3) the name of a function (e.g. 'd_fun') that takes
%            one argument, i, and returns a row vector containng the 
%            distances from all N points to point i. 
%
%    n_fcn = neighborhood function ('epsilon' or 'k') 
%    n_size = neighborhood size (value for epsilon or k) 
%
%    options = optional structure of options:
%      options.dims = (row) vector of embedding dimensionalities to use
%                        (1:10 = default)
%      options.comp = which connected component to embed, if more than one. 
%                        (1 = largest (default), 2 = second largest, ...)
%      options.display = plot residual variance and 2-D embedding?
%                        (1 = yes (default), 0 = no)
%      options.overlay = overlay graph on 2-D embedding?  
%                        (1 = yes (default), 0 = no)
%      options.verbose = display progress reports? 
%                        (1 = yes (default), 0 = no)
%      options.dijkstra = use dijkstra's algorithm for shortest paths with
%                         full N x N distance matrix. 
%                         (1 = yes (default), 0 = use Floyd; Floyd should
%                          be used only if you are unable to MEX dijkstra.cpp)
%      options.Kmax = maximum number of neighbors (used for sparse versions
%                        of epsilon; by default, estimated by random sample)
%      options.landmarks = (row) vector of landmark points to use in MDS. 
%                 (MDS finds the configuration that best approximates
%                  the distances from all points to the landmark points.
%                  The default landmark points are 1:N (i.e. all the points), 
%                  which is equivalent to classical MDS.  Good 
%                  results may often be obtained using a number of
%                  landmarks that is much smaller than N, but much
%                  larger than the data's intrinsic dimensionality.
%                  Note that this extension is experimental!  For
%                  discussion, see Steyvers, de Silva, and Tenenbaum
%                  (in preparation).)
%
% Output: 
%    Y = Y.coords is a cell array, with coordinates for d-dimensional embeddings
%         in Y.coords{d}.  Y.index contains the indices of the points embedded.
%    R = residual variances for embeddings in Y
%    E = edge matrix for neighborhood graph
%

%    BEGIN COPYRIGHT NOTICE
%
%    Isomap II code -- (c) 1998-2000 Josh Tenenbaum
%
%    This code is provided as is, with no guarantees except that 
%    bugs are almost surely present.  Published reports of research 
%    using this code (or a modified version) should cite the 
%    article that describes the algorithm: 
%
%      J. B. Tenenbaum, V. de Silva, J. C. Langford (2000).  A global
%      geometric framework for nonlinear dimensionality reduction.  
%      Science 290 (5500): 2319-2323, 22 December 2000.  
%
%    Comments and bug reports are welcome.  Email to jbt@psych.stanford.edu. 
%    I would also appreciate hearing about how you used this code, 
%    improvements that you have made to it, or translations into other
%    languages.    
%
%    You are free to modify, extend or distribute this code, as long 
%    as this copyright notice is included whole and unchanged.  
%
%    END COPYRIGHT NOTICE
%
% Modified by Ramon Casero <rcasero@gmail.com>, University of Oxford,  15
% Mar 2011 to fix some syntax problems.
%
% This file is distributed as a derivative work of a third-party function
% with project Gerardus.
%
% http://code.google.com/p/gerardus/
%
% Original code downloaded from the Isomap Homepage
%
% http://isomap.stanford.edu/


%%%%% Step 0: Initialization and Parameters %%%%%

if nargin < 3
     error('Too few input arguments'); 
elseif nargin < 4
     options = struct('dims',1:10,'overlay',1,'comp',1,'display',1,'dijkstra',1,'verbose',1); 
end

if ischar(D)
     mode = 3; 
     d_func = D; 
     N = length(feval(d_func,1)); 
elseif issparse(D) 
     mode = 2; 
     N = size(D,1); 
     if ~(N==size(D,2))
          error('D must be a square matrix'); 
     end; 
else 
     mode = 1; 
     N = size(D,1); 
     if ~(N==size(D,2))
         error('D must be a square matrix'); 
     end; 
end

if strcmp(n_fcn, 'k')
     K = n_size; 
     if ~(K==round(K))
         error('Number of neighbors for k method must be an integer');
     end
     if ((mode==2) && ~(min(sum(D'>0))>=K))
         error('Sparse D matrix must contain at least K nonzero entries in each row');
     end
elseif strcmp(n_fcn, 'epsilon')
     epsilon = n_size; 
     if isfield(options,'Kmax')
         K = options.Kmax; 
     elseif (mode==3)    %% estimate maximum equivalent K %% 
         tmp = zeros(10,N); 
         for i=1:10
             tmp(i,:) = feval(d_func,ceil(N*rand)); 
         end
         K = 2*max(sum(tmp'<epsilon));    % just to be safe
     end
else 
     error('Neighborhood function must be either epsilon or k'); 
end

if (mode == 3)
     INF = inf; 
else
     INF =  1000*max(max(D))*N;  %% effectively infinite distance
end

if ~isfield(options,'dims')
     options.dims = 1:10; 
end
if ~isfield(options,'overlay')
     options.overlay = 1; 
end
if ~isfield(options,'comp')
     options.comp = 1; 
end
if ~isfield(options,'display')
     options.display = 1; 
end
if ~isfield(options,'verbose')
     options.verbose = 1; 
end
if ~isfield(options,'landmarks')
     options.landmarks = 1:N; 
end
if ~isfield(options,'dijkstra')
     options.dijkstra = 1; 
end
dims = options.dims; 
comp = options.comp; 
overlay = options.overlay; 
displ = options.display; 
verbose = options.verbose; 
landmarks = options.landmarks; 
use_dijk = options.dijkstra;

Y.coords = cell(length(dims),1); 
R = zeros(1,length(dims)); 

%%%%% Step 1: Construct neighborhood graph %%%%%
disp('Constructing neighborhood graph...'); 

if ((mode == 1) && (use_dijk == 0))
     if strcmp(n_fcn, 'k')
         [tmp, ind] = sort(D); 
         tic; 
         for i=1:N
             D(i,ind((2+K):end,i)) = INF; 
             if ((verbose == 1) && (rem(i,50) == 0)) 
                 disp([' Iteration: ' num2str(i) '     Estimated time to completion: ' num2str((N-i)*toc/60/50) ' minutes']); tic; 
             end
         end
     elseif strcmp(n_fcn, 'epsilon')
         warning off    %% Next line causes an unnecessary warning, so turn it off
         D =  D./(D<=epsilon); 
         D = min(D,INF); 
         warning on
     end
     D = min(D,D');    %% Make sure distance matrix is symmetric
elseif ((mode == 1) && (use_dijk == 1))
     if n_fcn == 'k'
         [tmp, ind] = sort(D); 
         tic;
         for i=1:N
             D(i,ind((2+K):end,i)) = 0; 
             if ((verbose == 1) && (rem(i,50) == 0)) 
                 disp([' Iteration: ' num2str(i) '     Estimated time to completion: ' num2str((N-i)*toc/60/50) ' minutes']); tic; 
             end
         end
     elseif strcmp(n_fcn, 'epsilon')
         D =  D.*(D<=epsilon); 
     end
     D = sparse(D); 
     D = max(D,D');    %% Make sure distance matrix is symmetric
elseif (mode == 2)
     if n_fcn == 'k'
         Di = zeros(N*K,1);      Dj = zeros(N*K,1);       Ds = zeros(N*K,1); 
         counter = 0; 
         [a,b,c] = find(D); 
         tic; 
         for i=1:N
             l = find(a==i); 
             [g,f] = sort(c(l)); 
             Di(counter+(1:K)) = i; 
             Dj(counter+(1:K)) = b(l(f(1:K))); 
             Ds(counter+(1:K)) = g(1:K); 
             counter = counter+K; 
             if ((verbose == 1) && (rem(i,50) == 0)) 
                  disp([' Iteration: ' num2str(i) '     Estimated time to completion: ' num2str((N-i)*toc/60/50) ' minutes']); tic; 
             end
         end
         D = sparse(Di(1:counter), Dj(1:counter), Ds(1:counter));
         clear Di Dj Ds counter; 
     elseif strcmp(n_fcn, 'epsilon')
         D =  D.*(D<=epsilon); 
     end
     D = max(D,D');    %% Make sure distance matrix is symmetric
elseif (mode == 3)
     Di = zeros(N*(K+1),1);      Dj = zeros(N*(K+1),1);       Ds = zeros(N*(K+1),1); 
     counter = 0; 
     tic; 
     for i=1:N
         d = feval(d_func,i); 
         if n_fcn == 'k'
             [c,b] = sort(d); 
             Di(counter+(1:(K+1))) = i; 
             Dj(counter+(1:(K+1))) = b(1:(K+1)); 
             Ds(counter+(1:(K+1))) = c(1:(K+1)); 
             counter = counter+(K+1); 
         elseif strcmp(n_fcn, 'epsilon')
             [a,b,c] = find(d.*(d<=epsilon)); 
             l = length(a); 
             Di(counter+(1:l)) = i; 
             Dj(counter+(1:l)) = b; 
             Ds(counter+(1:l)) = c; 
             counter = counter+l; 
         end
         if ((verbose == 1) && (rem(i,50) == 0)) 
              disp([' Iteration: ' num2str(i) '     Estimated time to completion: ' num2str((N-i)*toc/60/50) ' minutes']); tic; 
         end
     end
     D = sparse(Di(1:counter), Dj(1:counter), Ds(1:counter));
     clear Di Dj Ds counter; 
     D = max(D,D');    %% Make sure distance matrix is symmetric
end    

if (overlay == 1)
     if ((mode == 1) && (use_dijk == 0))
         E = int8(1-(D==INF));  %%  Edge information for subsequent graph overlay
     else
         [a,b,c] = find(D); 
         E = sparse(a,b,ones(size(a))); 
     end
end

%%%%% Step 2: Compute shortest paths %%%%%
disp('Computing shortest paths...'); 

if ((mode==1) && (use_dijk == 0))
     tic; 
     for k=1:N
         D = min(D,repmat(D(:,k),[1 N])+repmat(D(k,:),[N 1])); 
         if ((verbose == 1) && (rem(k,20) == 0)) 
              disp([' Iteration: ' num2str(k) '     Estimated time to completion: ' num2str((N-i)*toc/i/60) ' minutes']); 
         end
     end
else
     D = dijkstra(D, landmarks);
end

%%%%% Step 3: Construct low-dimensional embeddings (Classical MDS) %%%%%
disp('Constructing low-dimensional embeddings (Classical MDS)...'); 

%%%%% Remove outliers from graph %%%%%
disp('  Checking for outliers...'); 

if ((mode == 1) && (use_dijk == 0))
     [tmp, firsts] = min(D==INF);     %% first point each point connects to
else
     [tmp, firsts] = min(D==inf);     %% first point each point connects to
end
[comps, I, J] = unique(firsts);    %% first point in each connected component
n_comps = length(comps);           %% number of connected components
size_comps = sum((repmat(firsts,n_comps,1)==((1:n_comps)'*ones(1,N)))'); 
                                   %% size of each connected component
[tmp, comp_order] = sort(size_comps);  %% sort connected components by size
comps = comps(comp_order(end:-1:1));    
size_comps = size_comps(comp_order(end:-1:1)); 
if (comp>n_comps)                
     comp=1;                              %% default: use largest component
end
Y.index = find(firsts==comps(comp)); %% list of points in relevant component
Y.index = setdiff(Y.index,find(isinf(min(D)))); %% prune points that don't connect
                                                %% to any landmarks
N = length(Y.index); 
[tmp, landmarks, land_ind] = intersect(landmarks,Y.index); 
                                       %% list of landmarks in component
nl = length(landmarks); 
D = full(D(landmarks,Y.index))'; 
disp(['    Number of connected components in graph: ' num2str(n_comps)]); 
disp(['    Embedding component ' num2str(comp) ' with ' num2str(length(Y.index)) ' points.']); 

dims = unique(min(dims,nl-1));    %% don't embed in more dimensions than landmarks-1
if (nl==N)
     opt.disp = 0; 
     [vec, val] = eigs(-.5*(D.^2 - sum(D.^2)'*ones(1,N)/N - ones(N,1)*sum(D.^2)/N + sum(sum(D.^2))/(N^2)), max(dims), 'LR', opt); 
else
     subB = -.5*(D.^2 - sum(D'.^2)'*ones(1,nl)/nl - ones(N,1)*sum(D.^2)/N+sum(sum(D.^2))/(N*nl));
     opt.disp = 0; 
     [alpha,beta] = eigs(subB'*subB, max(dims), 'LR', opt); 
     val = beta.^(1/2); 
     vec = subB*alpha*inv(val); 
     clear subB alpha beta; 
end
h = real(diag(val)); 
[foo,sorth] = sort(h);  sorth = sorth(end:-1:1); 
val = real(diag(val(sorth,sorth))); 
vec = vec(:,sorth); 

D = reshape(D,N*nl,1); 
for di = 1:length(dims)
     Y.coords{di} = real(vec(:,1:dims(di)).*(ones(N,1)*sqrt(val(1:dims(di)))'))'; 
     r2 = 1-corrcoef(reshape(real(L2_distance(Y.coords{di}, Y.coords{di}(:,land_ind))),N*nl,1),D).^2; 
     R(di) = r2(2,1); 
     if (verbose == 1)
         disp(['  Isomap on ' num2str(N) ' points with dimensionality ' num2str(dims(di)) '  --> residual variance = ' num2str(R(di))]); 
     end
end

clear D; 

%%%%%%%%%%%%%%%%%% Graphics %%%%%%%%%%%%%%%%%%

if (displ==1)
     %%%%% Plot fall-off of residual variance with dimensionality %%%%%
     figure;
     hold on
     plot(dims, R, 'bo'); 
     plot(dims, R, 'b-'); 
     hold off
     ylabel('Residual variance'); 
     xlabel('Isomap dimensionality'); 

     %%%%% Plot two-dimensional configuration %%%%%
     twod = find(dims==2); 
     if ~isempty(twod)
         figure;
         hold on;
         plot(Y.coords{twod}(1,:), Y.coords{twod}(2,:), 'ro'); 
         if (overlay == 1)
             gplot(E(Y.index, Y.index), [Y.coords{twod}(1,:); Y.coords{twod}(2,:)]'); 
             title('Two-dimensional Isomap embedding (with neighborhood graph).'); 
         else
             title('Two-dimensional Isomap.'); 
         end
         hold off;
     end
end

return;







