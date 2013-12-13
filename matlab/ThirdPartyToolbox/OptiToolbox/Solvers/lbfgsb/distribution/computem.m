function m = computem (L, w, W, K, Phi)

  % Get the number of documents.
  D = length(L);  

  % Initialize the counts "m".
  m = zeros(W,K);
  
  % Repeat for every document, then for every word in the document.
  is = 0;
  for d = 1:D
    is = is(end) + (1:L(d));
    for i = is
      j      = w(i);
      m(j,:) = m(j,:) + Phi(:,i)';
    end
  end

