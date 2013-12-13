function [xi, gamma, Phi] = mfldainit (W, K, w)
  
  % Get the number of documents.
  D = length(w);
  
  % Get the size of the corpus.
  M = 0;
  for d = 1:D
    M = M + length(w{d});
  end
  
  % Generate a random starting point for xi.
  xi = rand(W,K);
  
  % Generate a random starting point for gamma.
  gamma = rand(K,D);
  
  % Generate a random starting point for phi.
  Phi = rand(K,M);
  Phi = Phi ./ repmat(sum(Phi),K,1);
