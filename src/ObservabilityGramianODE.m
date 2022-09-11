function WdotVec = ObservabilityGramianODE(Wvec,A,C)
% Observability gramian ode equation for ode solver

% convert W to matrix
W = reshape(Wvec,size(A));

% compute Wdot
Wdot = -A'*W - W*A - C'*C;

% conver Pdot to vector
WdotVec = Wdot(:);
    
end