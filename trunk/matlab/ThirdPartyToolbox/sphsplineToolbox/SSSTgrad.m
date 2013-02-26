function y = SSSTgrad (x, p)
% SSSTgrad  Spherical Surface Spline in Tension gradient
%
%	y = SSSTgrad (x, p)
%
%   x is cos(theta) in range -1 <= x <= 1
%   p is tension (p >= 0)
%
% Returns the gradient of the Green's' function for a spherical
% surface spline in tension, following Wessel and Becker [2008].
% If p == 0 or not given then we return Parker's [1994] minimum
% curvature solution instead.

% $Id: SSSTgrad.m,v 1.1.1.1 2008/05/09 21:34:52 myself Exp $
% P. Wessel, SOEST, U of Hawaii, April 2008 (pwessel@hawaii.edu)

v = tension2nu (p);
y = zeros(size(x));
k = find (abs(x) < 1);
if (v == 0)     % Just do Parker's dilog solution
    y(k) = 0.5 * log(0.5 + 0.5 * x(k)) .* sqrt ((1.0 + x(k)) ./ (1.0 - x(k)));
    return
end
v1 = v + 1;
b = pi * v1 / sin(v*pi);
s = sqrt (1 - x(k).^2);
y(k) = -s.*(b*(x(k).*Pv(-x(k),v) + Pv(-x(k),v1))./(1-x(k).^2) + (1-x(k)).^(-1));
% s = sqrt (1 - x(k).^2);
% y(k) = -(b*(x(k).*Pv(-x(k),v) + Pv(-x(k),v1))./s + (1-x(k)).^(-1));
y = real(y);

% Sub-functions used by SSSTgrad
function P = Pv (x, v)

P = zeros(size(x));
for i = 1:length(x)
    if (x(i) == -1)
        p = inf;
    else
        [p q k] = PvQv (x(i), v);
    end
    P(i) = p;
end

function [Pv Qv iter] = PvQv (x, v)
% Based on recipe in "An Atlas of Functions" by
% Spanier and Oldham, 1987
iter = 0;
if (x == -1)
    Pv = -inf;
    Qv = -inf;
    return
end
if (x == +1)
    Pv = 1;
    Qv = inf;
    return
end
a = 1;
R = 1;
K = 4 * sqrt (abs(v - v^2));
if (abs(1+v) + floor (1+v)) == 0
	a = 1.0e99;
	v = -1 - v;
end
s = sin (0.5*pi*v);
c = cos (0.5*pi*v);
w = (0.5 + v)^2;
while v <= 6.0
	v = v + 2;
	R = R * (v - 1)/v;
end
X = 1.0 / (4 + 4*v);
g = 1 + 5*X*(1 - 3*X*(0.35+6.1*X));
R = R*(1 - X*(1 - g*X/2))/sqrt (8*X);
g = 2*x;
u = g;
f = 1;
t = 1;
k = 0.5;
X = 1 + (1e8/(1 - x.^2));

t = t .* x.^2 * (k^2 - w) / ((k + 1)^2 - 0.25);
k = k + 1;
f = f + t;
u = u .* x.^2 * (k^2 - w) / ((k + 1)^2 - 0.25);
k = k + 1;
g = g + u;
while (k < K || abs (X*t) > abs(f))
        iter = iter + 1;
	t = t .* x.^2 * (k^2 - w) / ((k + 1)^2 - 0.25);
	k = k + 1;
	f = f + t;
	u = u .* x.^2 * (k^2 - w) / ((k + 1)^2 - 0.25);
	k = k + 1;
	g = g + u;
end
f = f + (x.^2.*t ./ (1 - x.^2));
g = g + (x.^2.*u ./ (1 - x.^2));
Pv = ((s*g*R) + (c*f/R))/sqrt(pi);
Qv = a*sqrt(pi)*((c*g*R) - (s*f/R))/2;

function [f] = psi(z)
%Psi     Psi (or Digamma) function valid in the entire complex plane.
%
%                 d
%        Psi(z) = --log(Gamma(z))
%                 dz
%
%usage: [f] = psi(z)
%
%        Z may be complex and of any size.
%
%        This program uses the analytical derivative of the
%        Log of an excellent Lanczos series approximation
%        for the Gamma function.
%        
%References: C. Lanczos, SIAM JNA  1, 1964. pp. 86-96
%            Y. Luke, "The Special ... approximations", 1969 pp. 29-31
%            Y. Luke, "Algorithms ... functions", 1977
%            J. Spouge,  SIAM JNA 31, 1994. pp. 931
%            W. Press,  "Numerical Recipes"
%            S. Chang, "Computation of special functions", 1996
%
siz = size(z);
z=z(:);
zz=z;

f = zeros(size(z));

%reflection point
p=find(real(z)<0.5);
if ~isempty(p)
   z(p)=1-z(p);
end

%Lanczos approximation for the complex plane
 
g=607/128; % best results when 4<=g<=5
 
c = [  0.99999999999999709182;
      57.156235665862923517;
     -59.597960355475491248;
      14.136097974741747174;
      -0.49191381609762019978;
        .33994649984811888699e-4;
        .46523628927048575665e-4;
       -.98374475304879564677e-4;
        .15808870322491248884e-3;
       -.21026444172410488319e-3;
        .21743961811521264320e-3;
       -.16431810653676389022e-3;
        .84418223983852743293e-4;
       -.26190838401581408670e-4;
        .36899182659531622704e-5];


n=0;
d=0;
for k=size(c,1):-1:2
    dz=1./(z+k-2);
    dd=c(k).*dz;
    d=d+dd;
    n=n-dd.*dz;
end
d=d+c(1);
gg=z+g-0.5;
%log is accurate to about 13 digits...

f = log(gg) + (n./d - g./gg) ;

if ~isempty(p)
   f(p) = f(p)-pi*cot(pi*zz(p));
end

p=find(round(zz)==zz & real(zz)<=0 & imag(zz)==0);
if ~isempty(p)
   f(p) = Inf;
end

f=reshape(f,siz);

function y = dilog (x)
% DILOG   The dilogarithm
%
%   y = dilog (x)
%
% Compute dilog(x) (defined for x >= 0) by the method of Parker's
% Appendix A of his Geophysical Inverse Theory.  The function
% is needed for x in the range 0 <= x <= 1 when solving the
% spherical spline interpolation in section 2.07 of Parker.

y = zeros (size (x));

pisqon6 = pi * pi / 6.0;
k = find (x <= 0.0);
if (~isempty(k))
    y(k) = pisqon6;
end
k = find (x > 0.0 & x < 0.5);
if (~isempty(k))
	y(k) = -log (1.0 - x(k));
	ysq = y(k) .* y(k);
	z = y(k) .* (1.0 + y(k) .* (-0.25 + y(k) .* (0.027777777777213 + ...
        ysq .* (-2.7777776990e-04 + ysq .* (4.724071696e-06 + ...
        ysq .* (-9.1764954e-08 + 1.798670e-09 .* ysq))))));
	y(k) = pisqon6 - z + y(k) .* log (x(k));
end
k = find (x >= 0.5 & x < 2.0);
if (~isempty(k))
    y(k) = -log (x(k));
	ysq = y(k) .* y(k);
	z = y(k) .* (1.0 + y(k) .* (-0.25 + y(k) .* (0.027777777777213 + ...
        ysq .* (-2.7777776990e-04 + ysq .* (4.724071696e-06 + ...
        ysq .* (-9.1764954e-08 + 1.798670e-09 .* ysq))))));
    y(k) = z;
end
k = find (x >= 2.0);
if (~isempty(k))
    y(k) = log (x(k));
	ysq = y(k) .* y(k);
	z = y(k) .* (1.0 + y(k) .* (-0.25 + y(k) .* (0.027777777777213 + ...
        ysq .* (-2.7777776990e-04 + ysq .* (4.724071696e-06 + ...
        ysq .* (-9.1764954e-08 + 1.798670e-09 .* ysq))))));
        y(k) = -z - 0.5 * ysq;
end

function nu = tension2nu (p)
nu = (-1 + sqrt (1 - 4*p^2))/2;
