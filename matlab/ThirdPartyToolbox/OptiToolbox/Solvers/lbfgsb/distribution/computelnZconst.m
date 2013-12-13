function lnZconst = computelnZconst (nu, eta, D)

  K = length(nu);   % The number of topics.
  W = length(eta);  % The size of the vocabulary.

  sumnu    = sum(nu);
  sumeta   = sum(eta);
  lnZconst = K*gammaln(sumeta) - K*sum(gammaln(eta)) ...
             + D*gammaln(sumnu) - D*sum(gammaln(nu));
  