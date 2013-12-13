function [lnZ, xi, gamma, Phi] = mflda (nu, eta, w, xi, gamma, Phi)
  maxiter   = 1e4;
  m_lbfgs   = 20;
  tolerance = 0.1;
  lb        = 1e-4;
  
  % Make sure nu and eta are column vectors.
  nu  = nu(:);
  eta = eta(:);
  
  K = length(nu);   % The number of topics.
  W = length(eta);  % The size of the vocabulary.
  D = length(w);    % The number of documents.

  % Get the document lengths.
  L = zeros(1,D);
  for d = 1:D
    L(d) = length(w{d});
  end
  
  % Get the size of the corpus.
  w = [w{:}];
  M = length(w);
    
  % We precompute constant terms in the lower bound on the log-partition
  % function.
  lnZconst = computelnZconst(nu,eta,D);
  
  % Run the limited-memory BFGS algorithm.
  iterCallback = 'callbackMFLDA2';
  [xi gamma Phi] = lbfgsb({xi gamma Phi},...
                          {repmat(lb,W,K)  repmat(lb,K,D)  repmat(lb,K,M)},...
                         {repmat(inf,W,K) repmat(inf,K,D) repmat(inf,K,M)},...
                          'computeObjectiveMFLDA','computeGradientMFLDA',...
                          {nu eta L w lnZconst},'callbackMFLDA',...
                          'm',m_lbfgs,'factr',1e-12,'pgtol',...
                          tolerance,'maxiter',maxiter);

  % Normalize the topic proportions Phi.
  Phi = Phi ./ repmat(sum(Phi),K,1);  

  % Compute the logarithm of the normalizing constant.
  lnZ = computelnZ(L,w,nu,eta,xi,gamma,Phi,lnZconst);
  