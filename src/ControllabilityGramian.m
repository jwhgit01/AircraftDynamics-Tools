function [t,W] = ControllabilityGramian(A,B,tSpan)
%ControllabilityGramian

% make A,C functions of time
Afun = isa(A,'function_handle');
Bfun = isa(B,'function_handle');
if Afun
    At = A;
else
    At = @(t) A;
end
if Bfun
    Bt = B;
else
    Bt = @(t) B;
end

% if they are both constant, compute with Lyapunov eqation
if ~Afun && ~Bfun
    t = Inf;
    W = lyap(A,B*B');
    return
end

% size of system
n = size(At(0),1);

% gramian ODE
gramODE = @(t,Wvec) ControllabilityGramianODE(Wvec,At(t),Bt(t));

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