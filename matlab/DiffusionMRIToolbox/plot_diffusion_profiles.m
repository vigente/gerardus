function [ p, p2 ] = plot_diffusion_profiles( S, dash )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

x = 0:1E-5:3E-3;

if nargin == 1
    dash = 0;
end

if dash
    dd = '--';
else
    dd = '-';
end

hold on;


if isfield(S, 'strexp')
    u = S.strexp.D;
    beta = S.strexp.gamma;
    
    k = 0:5000;
    
    tww = u^(1/beta);
    l = zeros(size(x));
    for t_counter = 1:length(x)
    
        t = x(t_counter);
        
        % see wikipedia page on stretched exponential function
        %tosum = ((-1).^k)./factorial(k).*sin(pi*beta*k).*gamma(beta*k+1).*(t/tww).^(beta*k+1);
        tosum = ((-1).^k).*sin(pi*beta*k).*exp(gammaln(beta*k + 1) - gammaln(k + 1) + (beta*k+1)*log(t/tww));
        tosum(isnan(tosum)) = 0;

        l(t_counter) = -tww/(pi*t)*sum(tosum);
    end
    Y = l./x;
    Y(x < 0) = 0;
    Y(Y <= 0) = 0;
    Y(x > 2.5E-3) = 0; % to avoid numerical errors
    plot(x*1000, Y, dd, 'Color', [0.5 0.5 0.5], 'LineWidth', 2)

    
end

if isfield(S, 'DKI')
    D = S.DKI.D;
    K = S.DKI.K;
    
    sigma = sqrt(1/3 * D^2 * K);
    Y = normpdf(x, D, sigma);
    
    plot(x*1000, Y, [dd 'g'], 'LineWidth', 2);
    
end


if isfield(S, 'stat')
    mu = S.stat.mu;
    sigma = S.stat.sigma;
    if isfield(S.stat, 'Dmin')
        Dmin = S.stat.Dmin;
    else
        Dmin = 0;
    end
    if isfield(S.stat, 'Dmax')
        Dmax = S.stat.Dmax;
    else
        Dmax = inf;
    end
    
    Y = normpdf(x, mu, sigma) ...
      / (normcdf(Dmax, mu, sigma) - normcdf(Dmin, mu, sigma));
    Y(x < Dmin) = 0;
    Y(x > Dmax) = 0;
    
    plot(x*1000, Y, [dd 'k'], 'LineWidth', 2);
    
end


if isfield(S, 'gamma')
%     D = S.gamma.D;
%     sigma = S.gamma.sigma;
%     kappa = D^2 / sigma^2;
%     theta = sigma^2 / D;
    
    kappa = S.gamma.kappa;
    theta = S.gamma.theta;

    Ygamma = x.^(kappa - 1) .* exp(- x / theta) ./ (gamma(kappa) * theta^kappa);
    Ygamma(x < 0) = 0;
    
    plot(x*1000, Ygamma, dd, 'Color', [1 0.69 0.39], 'LineWidth', 2);
    
end


if 0%isfield(S, 'doubletrunc')
    
    mu = S.doubletrunc.mu;
    sigma = S.doubletrunc.sigma;
    D_upper = S.doubletrunc.D_upper;

    Y = normpdf(x, mu, sigma) / (normcdf(D_upper, mu, sigma) - normcdf(0, mu, sigma));
    Y(x < 0) = 0;
    Y(x > D_upper) = 0;

    plot(x*1000, Y, 'k--', 'LineWidth', 2);
    
end


if isfield(S, 'beta')
    a = S.beta.alpha;
    b = S.beta.beta;
    D0 = S.beta.D0;

    P = 1/D0 * gamma(a + b) / (gamma(a) * gamma(b)) *  (x/D0).^(a - 1) .* (1 - x / D0).^(b - 1);
    P(x < 0) = 0;
    P(x > D0) = 0;
    
    p = plot(x*1000, real(P), [dd 'r'], 'LineWidth', 2);
end



if isfield(S, 'monoexp')
    D = S.monoexp;
    plot([D D]*1000, [0 1500], dd, 'LineWidth', 2);
    
end


if isfield(S, 'biexp')
    Dfast = S.biexp.Dfast;
    Dslow = S.biexp.Dslow;
    f = S.biexp.f;
    dd = '-';
    p2 = plot([Dfast Dfast NaN Dslow Dslow]*1000, [0 1500*f NaN 0 1500*(1-f)], [dd 'c'], 'LineWidth', 2);

end


%legend('Mono-exponential', 'Bi-exponential', 'DKI (Gaussian)', 'Truncated Gaussian', 'Gamma', 'Beta')

set(gca, 'FontSize', 30)
set(gca, 'LineWidth', 2)
set(gca, 'XLim', [0 3])

