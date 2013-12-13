% GAMMARND(A) produces a single random deviate from the Gamma
% distribution with mean A and variance A.

function x = gammarnd (a)

  reject = true;
  while reject
    y0     = log(a)-1/sqrt(a);
    c      = a - exp(y0);
    y      = log(rand).*sign(rand-0.5)/c + log(a);
    f      = a*y-exp(y) - (a*y0 - exp(y0));
    g      = c*(abs((y0-log(a))) - abs(y-log(a)));
    reject = (log(rand) + g) > f;
  end

  x = exp(y);
