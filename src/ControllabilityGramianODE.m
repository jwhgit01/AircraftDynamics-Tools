function WdotVec = ControllabilityGramianODE(Wvec,A,B)
% Controllability gramian ode equation for ode solver

% convert W to matrix
W = reshape(Wvec,size(A));

% compute Wdot
Wdot = -A*W - W*A' - B*B';

% conver Pdot to vector
WdotVec = Wdot(:);
    
end