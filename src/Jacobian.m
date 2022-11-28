function J = Jacobian(f,x,Method,epsilon)
% Calculate Jacobian of function f at given x

% tolerance
if nargin < 4
    epsilon = 1e-6; 
end
epsilon_inv = 1/epsilon;

% evaluate f at x
f0 = feval(f,x); 

% dimensions
n = length(x);
m = length(f0);

%
% Single Step Finite Difference
%
if strcmp(Method,'FiniteDifference')
    J = zeros(m,n);
    for ii = 1:n
        xplus = x;
        xplus(ii) =  x(ii) + epsilon;
        J(:, ii) = (feval(f, xplus)-f0)*epsilon_inv;
    end
%
% Complex-Step Derivative Approximation (CSDA)
%
elseif strcmp(Method,'CSDA')
    J = zeros(m,n);
    for ii = 1:n
        xplus = x;
        xplus(ii) =  x(ii) + 1i*epsilon;
        J(:, ii) = imag(feval(f,xplus)).*epsilon_inv;
    end
%
% Central Finite Difference Method (CFDM)
%
elseif strcmp(Method,'CFDM')
    J = zeros(m,n);
    for ii = 1:n
        xplus = x;
        xminus = x;
        xplus(ii) =  x(ii) + epsilon;
        xminus(ii) =  x(ii) - epsilon;
        J(:, ii) = 0.5*epsilon_inv*(feval(f,xplus)-feval(f, xminus));
    end
%
% Symbolic computation
%
elseif strcmp(Method,'Symbolic')
    syms xs [n 1]
    fs = feval(f,xs);
    fx = jacobian(fs,xs);
    J = double(subs(fx,xs,x));
else
    error([Method ' method not recognized!']);
end
