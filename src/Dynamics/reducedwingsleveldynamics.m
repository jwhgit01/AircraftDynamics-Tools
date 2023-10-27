function f = reducedwingsleveldynamics(z,params)

% Parse input
theta = z(1,1);
u = z(2,1);
w = z(3,1);
de = z(4,1);
Omega = params.Omega;

% Aircraft and enviromental Constants
m = params.m;
I = params.I;
g = params.g;
rho = params.rho;

% Rotation matrix
psi = 0;
phi = 0;
R1 = [1,0,0;0,cos(phi),-sin(phi);0,sin(phi),cos(phi)];
R2 = [cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)];
R3 = [cos(psi),-sin(psi),0;sin(psi),cos(psi),0;0,0,1];
R_IB = R3*R2*R1;

% Aerodynamic forces and moments
vb = [u;0;w];
omega = zeros(3,1);
delta = [0;de;0;Omega];
Parameters = params.ParameterEstimates;
Constants.rho = rho;
Constants.g = g;
[F,M] = params.AerodynamicModel.Model([],vb,omega,delta,Parameters,Constants);

% Altitude kinematics
qdot = R_IB*vb;
zdot = qdot(3,1);
            
% Translational dynamics
e3 = [0;0;1];
vbdot = cross(vb,omega) + g*R_IB'*e3 + F/m;
omegadot = I\(cross(I*omega,omega) + M);

% Function for which we want to find zeros
f = [zdot;vbdot;omegadot];

end