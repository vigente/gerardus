function r = dirichletrnd (a)
  
  if nargin < 2
    uselightspeed = false;
  end

  % Take a sample from a Dirichlet distribution
  n = length(a);
  r = zeros(1,n);
  for i = 1:n
    r(i) = gammarnd(a(i));
  end
  r = r / sum(r);
