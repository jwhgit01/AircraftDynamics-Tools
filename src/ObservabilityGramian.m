function [t,W] = ObservabilityGramian(A,C,tSpan)
%ObservabilityGramian

% make A,C functions of time
Afun = isa(A,'function_handle');
Cfun = isa(C,'function_handle');
if Afun
    At = A;
else
    At = @(t) A;
end
if Cfun
    Ct = C;
else
    Ct = @(t) C;
end

% if they are both constant, compute with Lyapunov eqation
if ~Afun && ~Cfun
    t = Inf;
    W = lyap(A',C'*C);
    return
end

% size of system
n = size(At(0),1);

% gramian ODE
gramODE = @(t,Wvec) ObservabilityGramianODE(Wvec,At(t),Ct(t));

% numerically solve backwards in time
options = odeset('AbsTol',1e-12,'RelTol',1e-6);
[tODE,WODE] = ode45(gramODE,fliplr(tSpan),zeros(n),options);

% number of time steps in solution
N = length(tODE);

% flip and reshape ode45 solution
t = flipud(tODE);
WODE = flipud(WODE);
W = zeros(N,n,n);
for ii = 1:N
    W(ii,:,:) = reshape(WODE(ii,:),[n,n]); 
end



end