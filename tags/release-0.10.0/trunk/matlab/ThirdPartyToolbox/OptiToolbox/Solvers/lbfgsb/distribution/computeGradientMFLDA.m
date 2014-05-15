function [dxi, dgamma, dPhi] = computeGradientMFLDA (xi, gamma, Phi, auxdata)
  [nu eta L w] = deal(auxdata{1:4});

  D = length(L);    % The number of documents.
  K = length(nu);   % The number of topics.
  W = length(eta);  % The size of the vocabulary.
  M = length(w);    % The size of the corpus.
  
  % Normalize the topic proportions Phi.
  c   = sum(Phi);
  Phi = Phi ./ repmat(c,K,1);  
  
  % Compute the partial first derivatives of the objective with respect
  % to the variational parameter xi.
  m   = computem(L,w,W,K,Phi);
  x   = repmat(eta,1,K) + m - xi;
  dxi = trigamma(xi).*x - repmat(trigamma(sum(xi)).*sum(x),W,1);

  % Compute the partial first derivatives of the objective with respect
  % to the variational parameter gamma. Repeat for each document.
  dgamma = zeros(K,D);
  is     = 0;
  for d = 1:D
    is          = is(end) + (1:L(d));
    phi         = Phi(:,is);
    g           = gamma(:,d);
    x           = nu + row_sum(phi) - g;
    dgamma(:,d) = trigamma(g).*x - repmat(trigamma(sum(g))*sum(x),K,1);
  end
  
  % Compute the partial first derivatives of the objective with respect
  % to the variational parameters phi. Repeat for each document, then for
  % each word in the document.
  dPhi = zeros(K,M);
  is   = 0;
  for d = 1:D
    is = is(end) + (1:L(d));
    g  = gamma(:,d);
    for i = is
      j         = w(i);
      dPhi(:,i) = digamma(xi(j,:))' - digamma(sum(xi))' + ...
                  digamma(g) - digamma(sum(g)) - log(Phi(:,i)) - 1;
    end
  end
  
  % We do need to fix one thing. Namely, we want the partial derivatives of
  % the objective with respect to the *unnormalized* variational parameters
  % phi, not their normalized counterparts.
  dPhi = (dPhi - repmat(sum(dPhi.*Phi),K,1)) ./ repmat(c,K,1);
  
  % What we have computed is actually the gradient of the variational
  % lower bound, but we need minus that since we're solving a
  % minimization problem.
  dxi    = -dxi;
  dgamma = -dgamma;
  dPhi   = -dPhi;

