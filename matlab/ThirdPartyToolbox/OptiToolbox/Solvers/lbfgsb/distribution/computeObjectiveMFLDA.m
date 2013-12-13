function f = computeObjectiveMFLDA (xi, gamma, Phi, auxdata)
  [nu eta L w lnZconst] = deal(auxdata{:});

  % Normalize the topic proportions Phi.
  K   = length(nu);
  Phi = Phi ./ repmat(sum(Phi),K,1);
  
  % Compute the logarithm of the normalizing constant.
  lnZ = computelnZ(L,w,nu,eta,xi,gamma,Phi,lnZconst);
  
  % The objective is the negative of lnZ since we are solving a
  % minimization problem, not a maximization problem.
  f = -lnZ;
  