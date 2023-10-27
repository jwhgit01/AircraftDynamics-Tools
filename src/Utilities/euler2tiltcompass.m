function [lambda,zeta] = euler2tiltcompass(Theta)
%euler2tiltcompass 

e1 = [1;0;0];
e3 = [0;0;1];

phi = Theta(1,1);
theta = Theta(2,1);
psi = Theta(3,1);

R1 = [1,0,0;0,cos(phi),-sin(phi);0,sin(phi),cos(phi)];
R2 = [cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)];
R3 = [cos(psi),-sin(psi),0;sin(psi),cos(psi),0;0,0,1];

R_BI = (R3*R2*R1).';

lambda = R_BI*e1;
zeta = R_BI*e3;

end