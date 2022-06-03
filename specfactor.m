function [h,q,SSdim] = specfactor(cov,theta,varargin)
%SPECFACTOR Decomposes covariance function using spectral factorization. 
% S(w) = H(iw)*q*H(-iw)
%
% [h,q,dim] = SPECFACTOR(cov, theta, N, p)
%   cov  : covariance function ('matern' or 'sqexp')
%   theta: vector of hyperparameters for covariance 'cov'
%          [sigma2, l]
%   N    : (opt) number of terms in sqexp approximation (must be even)
%   p    : (opt) smoothness parameter for matern functions

if strcmp(cov,'matern')
    % CURRENTLY ONLY WORKS FOR P=2
    
    % Matern covariance with p=2 (v=5/2)
    p = varargin{1};
    if nargout > 2
        SSdim = p+1;
    end
    sig2 = theta(1); 
    l = theta(2);
    
    %coefficients for input 
    lambda = sqrt(2*(p+1/2))/l;
    h0 = lambda^3;
    h1 = 3*lambda^2;
    h2 = 3*lambda;
    h = [h0, h1, h2];
    
    %spectral density 
    q = (2*sig2*sqrt(pi())*lambda^(2*p+1)*gamma(p+1))/gamma(p+1/2);
    return 
end

if strcmp(cov,'sqexp')
    % Squared Exponential (approx w/ N=6)
    sig2 = theta(1);
    l = theta(2);
    N = varargin{1};
    if mod(N,2) ~= 0
        error('N must be an even integer (N=%d)', N)
    end
    if nargout > 2
        SSdim = N;
    end
    
    k = 1/(2*l^2);
    
    Pd = zeros(1,2*N+1);
    for i=2*N+1:-2:1
        ii = (i-1)/2;
        Pd(i) = (-1)^ii*factorial(N)/factorial(ii)*(4*k)^(N-ii);
    end
    Pd = fliplr(Pd);
    rts = roots(Pd);
    % Pd2 should be equal to Pd (note:coefficients will be in reverse order)
    % Pd2 = double(coeffs((omega-rts(1))*(omega-rts(2))*(omega-rts(3))*(omega-rts(4))*(omega-rts(5))*(omega-rts(6))...
    %                     *(omega-rts(7))*(omega-rts(8))*(omega-rts(9))*(omega-rts(10))*(omega-rts(11))*(omega-rts(12))));

    %coefficients of P-, lowest order first
    syms omega
    h = double(coeffs((omega-rts(1))*(omega-rts(2))*(omega-rts(3))*(omega-rts(4))*(omega-rts(5))*(omega-rts(6))));
    h = h(1:N); 
    
    %spectral density
    q = sig2*factorial(N)*(4*k)^N*sqrt(pi/k);
    return 
end 

end

