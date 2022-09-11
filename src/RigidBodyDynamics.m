function dxdt = RigidBodyDynamics(x,F_body,M_body,mass,MomentOfInertia)
%RigidBodyDynamics
%
% The states of this system are
%   x = [x y z phi theta psi u v w p q r]
% with inputs of forces and moments

% enviromental constants
g = 9.81;

% aircraft states
Theta = x(4:6,1);
vb = x(7:9,1);
omega = x(10:12,1);

% euler angles
phi = Theta(1);
theta = Theta(2);
psi = Theta(3);

% rotation matrix functions
R1 = @(phi) [1,0,0;0,cos(phi),-sin(phi);0,sin(phi),cos(phi)];
R2 = @(theta) [cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)];
R3 = @(psi) [cos(psi),-sin(psi),0;sin(psi),cos(psi),0;0,0,1];

% rotation matrix from body to inertial frame and visa versa
R_IB = R3(psi)*R2(theta)*R1(phi);
R_BI = R_IB.';

% Thetadot = L_IB Omega
L_IB = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);...
        0, cos(phi),           -sin(phi)           ;...
        0, sin(phi)*sec(theta), cos(phi)*sec(theta)];

% translational kinematics
rvecdot = R_IB*vb;

% rotational kinematics
Thetadot = L_IB*omega;

% Linear and angular momentum
p_lm = mass*vb;
h = MomentOfInertia*omega;

% Dynamic equations of motion
e3 = [0;0;1];
vbdot = 1/mass*(cross(p_lm,omega) + mass*g*R_BI*e3 + F_body);
omegadot = MomentOfInertia\(cross(h,omega) + M_body);

% assemble output
dxdt = [rvecdot;Thetadot;vbdot;omegadot];

end