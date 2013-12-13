function x = osborne(p)
  n=33;

  for i=1:n
    t=10*(i-1);
    x(i)=p(1) + p(2)*exp(-p(4)*t) + p(3)*exp(-p(5)*t);
  end
