function [topics, w] = createsyntheticdata (nu, eta, K, W, L)
  
  % Get the number of documents to generate.
  D = length(L);
  
  % Generate the word proportions for each topic from the prior.
  topics = zeros(W,K);
  for k = 1:K
    topics(:,k) = dirichletrnd(eta);
  end
  
  % Generate the documents.
  w = cell(D,1);
  for d = 1:D
    
    % Generate the topic proportions from the prior.
    t = dirichletrnd(nu)';
    
    % Generate the topic associated with each word.
    k = sample(t,L(d));
    
    % Generate the word tokens.
    w{d} = int32(sample_vector(topics,k));
  end
  