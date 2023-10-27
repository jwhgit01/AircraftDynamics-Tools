function [xdot,A,D] = nonlindyn_windplusjerk(t,x,u,vtil,dervflag,params)
%nonlindyn
%
% Copyright (c) 2023 Jeremy W. Hopwood. All rights reserved.
%
% This is a function defines the aircraft dynamics in a frozen,
% uniform wind field in the form required by the Estimation-Tools
% toolbox. It is a model for the continuous-time, nonlinear model of
% a dynamical control system,
%
%                   dx/dt = f(t,x,u,vtil)                   (1)
%
% where vtil(t) is a continuous-time random process.
%
% Inputs:
%
%   t           The time at which dx/dt is evaluated.
%
%   x           The nx x 1 state vector at time t given by
%                   / Position               \
%                   | Attitude               |
%                   | Angular Velocity       |
%                   | Body Acceleration      |
%                   | Relative Body Velocity |
%                   \ Wind Velocity          /
%
%   u           The nu x 1 control vector at time t.
%
%   vtil        The nv x 1 process noise vector at time t.
%
%   dervflag    A flag that determines whether (dervflag = 1) or
%               not (dervflag = 0) the partial derivatives df/dx and
%               df/dvtil must be calculated. If dervflag = 0, then these
%               outputs will be empty arrays.
%
%   params      Any structure array that is used to pass constants through
%               to the model. The following fields are required:
%
%               g                  The positive gravitational
%                                  acceleration
%
%               rho                The constant air density
%
%               m                  The aircraft's mass in kg.
%
%               I                  The aircraft's moment of inertia.
%
%               AerodynamicModel   The handle to the desired
%                                  AerodynamicModel object
%
%               ParameterEstimates The handle to the desired
%                                  ParameterEstimates object or
%                                  numeric array of estimates
%
%               AttitudeParameterization The char array or string
%                                  that determines which attitude
%                                  paramterization is to be used.
%
%  
% Outputs:
%
%   xdot        The time derivative of x at time t from Eq. (1).
%
%   A           The Jacobian of f with respect to x. It is
%               evaluated and output only if dervflag = 1.
%               Otherwise, an empty array is output.
%
%   D           The Jacobian of f with respect to vtil. It is
%               evaluated and output only if dervflag = 1.
%               Otherwise, an empty array is output.
%

% Aircraft and enviromental Constants
m = params.m;
I = params.I;
g = params.g;

% Parse the state vector and compute attitude kinematics
Theta = x(4:6,1);
omega = x(7:9,1);
ab = x(10:12,1);
vbr = x(13:15,1);
w = x(16:18,1);

% Sources of process noise
if isempty(vtil)
    vtil = zeros(9,1);
end
vtil_Fadot = vtil(1:3,1); % time derivative of aerodynamic forces
vtil_Ma = vtil(4:6,1);    % aerodynamic moments
vtil_w = vtil(7:9,1);     % atmospheric turbulence

% Compute attitude kinematics
phi = Theta(1);
theta = Theta(2);
psi = Theta(3);
R1 = [1,0,0;0,cos(phi),-sin(phi);0,sin(phi),cos(phi)];
R2 = [cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)];
R3 = [cos(psi),-sin(psi),0;sin(psi),cos(psi),0;0,0,1];
R_IB = R3*R2*R1;
L_IB = [1, sin(phi)*tan(theta), cos(phi)*tan(theta);...
        0, cos(phi),           -sin(phi)           ;...
        0, sin(phi)*sec(theta), cos(phi)*sec(theta)];
attitudeDyn = L_IB*omega;

% Body velocity
vb = vbr + R_IB'*w;

% translational kinematics
NEDdot = R_IB*vbr + w;
            
% Dynamic equations of motion
vbrdot = cross(vbr,omega) + ab - R_IB'*vtil_w;
omegadot = I\(cross(I*omega,omega) + vtil_Ma);

% Translational jerk
e3 = [0;0;1];
abdot = cross(vb,omegadot) - g*cross(omega,R_IB'*e3) + cross(vb,I\vtil_Ma) + vtil_Fadot/m;

% Wind dynamics
wdot = vtil_w;

% Rigid body in wind state derivative
xdot = [NEDdot;attitudeDyn;omegadot;abdot;vbrdot;wdot];

% If the Jacobians are not needed, set to empty arrays and return
if dervflag == 0
    A = [];
    D = [];
    return
end

% Compute the Jacobains of the model-free aircraft dynamics
Ixx = I(1,1);
Iyy = I(2,2);
Izz = I(3,3);
Ixz = -I(1,3);
A = [0,0,0,x(14,1)*(sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1)))+x(15,1)*(cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1))),cos(x(6,1))*(x(15,1)*cos(x(4,1))*cos(x(5,1))-x(13,1)*sin(x(5,1))+x(14,1)*cos(x(5,1))*sin(x(4,1))),x(15,1)*(cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1)))-x(14,1)*(cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))-x(13,1)*cos(x(5,1))*sin(x(6,1)),0,0,0,0,0,0,cos(x(5,1))*cos(x(6,1)),cos(x(6,1))*sin(x(4,1))*sin(x(5,1))-cos(x(4,1))*sin(x(6,1)),sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1)),1,0,0;
    0,0,0,-x(14,1)*(cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1)))-x(15,1)*(cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1))),sin(x(6,1))*(x(15,1)*cos(x(4,1))*cos(x(5,1))-x(13,1)*sin(x(5,1))+x(14,1)*cos(x(5,1))*sin(x(4,1))),x(15,1)*(sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1)))-x(14,1)*(cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1)))+x(13,1)*cos(x(5,1))*cos(x(6,1)),0,0,0,0,0,0,cos(x(5,1))*sin(x(6,1)),cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)),cos(x(4,1))*sin(x(5,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1)),0,1,0;
    0,0,0,cos(x(5,1))*(x(14,1)*cos(x(4,1))-x(15,1)*sin(x(4,1))),-x(13,1)*cos(x(5,1))-x(15,1)*cos(x(4,1))*sin(x(5,1))-x(14,1)*sin(x(4,1))*sin(x(5,1)),0,0,0,0,0,0,0,-sin(x(5,1)),cos(x(5,1))*sin(x(4,1)),cos(x(4,1))*cos(x(5,1)),0,0,1;
    0,0,0,(sin(x(5,1))*(x(8,1)*cos(x(4,1))-x(9,1)*sin(x(4,1))))/cos(x(5,1)),(x(9,1)*cos(x(4,1))+x(8,1)*sin(x(4,1)))/cos(x(5,1))^2,0,1,sin(x(4,1))*tan(x(5,1)),cos(x(4,1))*tan(x(5,1)),0,0,0,0,0,0,0,0,0;
    0,0,0,-x(9,1)*cos(x(4,1))-x(8,1)*sin(x(4,1)),0,0,0,cos(x(4,1)),-sin(x(4,1)),0,0,0,0,0,0,0,0,0;
    0,0,0,(x(8,1)*cos(x(4,1))-x(9,1)*sin(x(4,1)))/cos(x(5,1)),(sin(x(5,1))*(x(9,1)*cos(x(4,1))+x(8,1)*sin(x(4,1))))/cos(x(5,1))^2,0,0,sin(x(4,1))/cos(x(5,1)),cos(x(4,1))/cos(x(5,1)),0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,-(Ixz*x(8,1)*(Ixx-Iyy+Izz))/(Ixz^2-Ixx*Izz),(Ixz^2*x(9,1)+Izz^2*x(9,1)-Ixx*Ixz*x(7,1)+Ixz*Iyy*x(7,1)-Ixz*Izz*x(7,1)-Iyy*Izz*x(9,1))/(Ixz^2-Ixx*Izz),(x(8,1)*(Ixz^2+Izz^2-Iyy*Izz))/(Ixz^2-Ixx*Izz),0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,-(Ixx*x(9,1)+2*Ixz*x(7,1)-Izz*x(9,1))/Iyy,0,(2*Ixz*x(9,1)-Ixx*x(7,1)+Izz*x(7,1))/Iyy,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,-(x(8,1)*(Ixx^2-Iyy*Ixx+Ixz^2))/(Ixz^2-Ixx*Izz),-(Ixx^2*x(7,1)+Ixz^2*x(7,1)-Ixx*Ixz*x(9,1)-Ixx*Iyy*x(7,1)+Ixz*Iyy*x(9,1)-Ixz*Izz*x(9,1))/(Ixz^2-Ixx*Izz),(Ixz*x(8,1)*(Ixx-Iyy+Izz))/(Ixz^2-Ixx*Izz),0,0,0,0,0,0,0,0,0;
    0,0,0,g*(x(9,1)*cos(x(4,1))*cos(x(5,1))+x(8,1)*cos(x(5,1))*sin(x(4,1)))-((Ixx*vtil(6,1)+Ixz*vtil(4,1))*(x(16,1)*(sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1)))-x(17,1)*(cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(18,1)*cos(x(4,1))*cos(x(5,1))))/(Ixz^2-Ixx*Izz)+(vtil(5,1)*(x(17,1)*(cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))-x(16,1)*(cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1)))+x(18,1)*cos(x(5,1))*sin(x(4,1))))/Iyy-((x(16,1)*(sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1)))-x(17,1)*(cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(18,1)*cos(x(4,1))*cos(x(5,1)))*(Ixx*vtil(6,1)+Ixz*vtil(4,1)+Ixx^2*x(7,1)*x(8,1)+Ixz^2*x(7,1)*x(8,1)-Ixx*Ixz*x(8,1)*x(9,1)-Ixx*Iyy*x(7,1)*x(8,1)+Ixz*Iyy*x(8,1)*x(9,1)-Ixz*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz)+((x(17,1)*(cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))-x(16,1)*(cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1)))+x(18,1)*cos(x(5,1))*sin(x(4,1)))*(vtil(5,1)-Ixz*x(7,1)^2+Ixz*x(9,1)^2-Ixx*x(7,1)*x(9,1)+Izz*x(7,1)*x(9,1)))/Iyy,g*sin(x(5,1))*(x(8,1)*cos(x(4,1))-x(9,1)*sin(x(4,1)))-(sin(x(4,1))*(x(16,1)*cos(x(5,1))*cos(x(6,1))-x(18,1)*sin(x(5,1))+x(17,1)*cos(x(5,1))*sin(x(6,1)))*(Ixx*vtil(6,1)+Ixz*vtil(4,1)+Ixx^2*x(7,1)*x(8,1)+Ixz^2*x(7,1)*x(8,1)-Ixx*Ixz*x(8,1)*x(9,1)-Ixx*Iyy*x(7,1)*x(8,1)+Ixz*Iyy*x(8,1)*x(9,1)-Ixz*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz)-(cos(x(4,1))*(x(16,1)*cos(x(5,1))*cos(x(6,1))-x(18,1)*sin(x(5,1))+x(17,1)*cos(x(5,1))*sin(x(6,1)))*(vtil(5,1)-Ixz*x(7,1)^2+Ixz*x(9,1)^2-Ixx*x(7,1)*x(9,1)+Izz*x(7,1)*x(9,1)))/Iyy-(sin(x(4,1))*(Ixx*vtil(6,1)+Ixz*vtil(4,1))*(x(16,1)*cos(x(5,1))*cos(x(6,1))-x(18,1)*sin(x(5,1))+x(17,1)*cos(x(5,1))*sin(x(6,1))))/(Ixz^2-Ixx*Izz)-(vtil(5,1)*cos(x(4,1))*(x(16,1)*cos(x(5,1))*cos(x(6,1))-x(18,1)*sin(x(5,1))+x(17,1)*cos(x(5,1))*sin(x(6,1))))/Iyy,((x(16,1)*(cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(17,1)*(cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1))))*(Ixx*vtil(6,1)+Ixz*vtil(4,1)+Ixx^2*x(7,1)*x(8,1)+Ixz^2*x(7,1)*x(8,1)-Ixx*Ixz*x(8,1)*x(9,1)-Ixx*Iyy*x(7,1)*x(8,1)+Ixz*Iyy*x(8,1)*x(9,1)-Ixz*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz)-((x(16,1)*(cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(17,1)*(sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1))))*(vtil(5,1)-Ixz*x(7,1)^2+Ixz*x(9,1)^2-Ixx*x(7,1)*x(9,1)+Izz*x(7,1)*x(9,1)))/Iyy+((x(16,1)*(cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(17,1)*(cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1))))*(Ixx*vtil(6,1)+Ixz*vtil(4,1)))/(Ixz^2-Ixx*Izz)-(vtil(5,1)*(x(16,1)*(cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(17,1)*(sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1)))))/Iyy,((Ixx*x(9,1)+2*Ixz*x(7,1)-Izz*x(9,1))*(x(15,1)+x(16,1)*(sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1)))-x(17,1)*(cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(18,1)*cos(x(4,1))*cos(x(5,1))))/Iyy-((x(8,1)*Ixx^2-Iyy*x(8,1)*Ixx+x(8,1)*Ixz^2)*(x(14,1)-x(16,1)*(cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1)))+x(17,1)*(cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(18,1)*cos(x(5,1))*sin(x(4,1))))/(Ixz^2-Ixx*Izz),-((x(14,1)-x(16,1)*(cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1)))+x(17,1)*(cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(18,1)*cos(x(5,1))*sin(x(4,1)))*(Ixx^2*x(7,1)+Ixz^2*x(7,1)-Ixx*Ixz*x(9,1)-Ixx*Iyy*x(7,1)+Ixz*Iyy*x(9,1)-Ixz*Izz*x(9,1)))/(Ixz^2-Ixx*Izz)-g*cos(x(4,1))*cos(x(5,1)),g*cos(x(5,1))*sin(x(4,1))-((2*Ixz*x(9,1)-Ixx*x(7,1)+Izz*x(7,1))*(x(15,1)+x(16,1)*(sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1)))-x(17,1)*(cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(18,1)*cos(x(4,1))*cos(x(5,1))))/Iyy+((Ixx*Ixz*x(8,1)-Ixz*Iyy*x(8,1)+Ixz*Izz*x(8,1))*(x(14,1)-x(16,1)*(cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1)))+x(17,1)*(cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(18,1)*cos(x(5,1))*sin(x(4,1))))/(Ixz^2-Ixx*Izz),0,0,0,0,-(2*Ixx*vtil(6,1)+2*Ixz*vtil(4,1)+Ixx^2*x(7,1)*x(8,1)+Ixz^2*x(7,1)*x(8,1)-Ixx*Ixz*x(8,1)*x(9,1)-Ixx*Iyy*x(7,1)*x(8,1)+Ixz*Iyy*x(8,1)*x(9,1)-Ixz*Izz*x(8,1)*x(9,1))/(Ixz^2-Ixx*Izz),-(2*vtil(5,1)-Ixz*x(7,1)^2+Ixz*x(9,1)^2-Ixx*x(7,1)*x(9,1)+Izz*x(7,1)*x(9,1))/Iyy,((cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1)))*(Ixx*vtil(6,1)+Ixz*vtil(4,1)+Ixx^2*x(7,1)*x(8,1)+Ixz^2*x(7,1)*x(8,1)-Ixx*Ixz*x(8,1)*x(9,1)-Ixx*Iyy*x(7,1)*x(8,1)+Ixz*Iyy*x(8,1)*x(9,1)-Ixz*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz)-(vtil(5,1)*(sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1))))/Iyy-((sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1)))*(vtil(5,1)-Ixz*x(7,1)^2+Ixz*x(9,1)^2-Ixx*x(7,1)*x(9,1)+Izz*x(7,1)*x(9,1)))/Iyy+((cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1)))*(Ixx*vtil(6,1)+Ixz*vtil(4,1)))/(Ixz^2-Ixx*Izz),(vtil(5,1)*(cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1))))/Iyy-((cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))*(Ixx*vtil(6,1)+Ixz*vtil(4,1)+Ixx^2*x(7,1)*x(8,1)+Ixz^2*x(7,1)*x(8,1)-Ixx*Ixz*x(8,1)*x(9,1)-Ixx*Iyy*x(7,1)*x(8,1)+Ixz*Iyy*x(8,1)*x(9,1)-Ixz*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz)+((cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1)))*(vtil(5,1)-Ixz*x(7,1)^2+Ixz*x(9,1)^2-Ixx*x(7,1)*x(9,1)+Izz*x(7,1)*x(9,1)))/Iyy-((cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))*(Ixx*vtil(6,1)+Ixz*vtil(4,1)))/(Ixz^2-Ixx*Izz),-(cos(x(5,1))*sin(x(4,1))*(Ixx*vtil(6,1)+Ixz*vtil(4,1)))/(Ixz^2-Ixx*Izz)-(vtil(5,1)*cos(x(4,1))*cos(x(5,1)))/Iyy-(cos(x(5,1))*sin(x(4,1))*(Ixx*vtil(6,1)+Ixz*vtil(4,1)+Ixx^2*x(7,1)*x(8,1)+Ixz^2*x(7,1)*x(8,1)-Ixx*Ixz*x(8,1)*x(9,1)-Ixx*Iyy*x(7,1)*x(8,1)+Ixz*Iyy*x(8,1)*x(9,1)-Ixz*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz)-(cos(x(4,1))*cos(x(5,1))*(vtil(5,1)-Ixz*x(7,1)^2+Ixz*x(9,1)^2-Ixx*x(7,1)*x(9,1)+Izz*x(7,1)*x(9,1)))/Iyy;
    0,0,0,((Ixz*vtil(6,1)+Izz*vtil(4,1))*(x(17,1)*(cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))-x(16,1)*(cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1)))+x(18,1)*cos(x(5,1))*sin(x(4,1))))/(Ixz^2-Ixx*Izz)+((x(17,1)*(cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))-x(16,1)*(cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1)))+x(18,1)*cos(x(5,1))*sin(x(4,1)))*(Ixz*vtil(6,1)+Izz*vtil(4,1)-Ixz^2*x(8,1)*x(9,1)-Izz^2*x(8,1)*x(9,1)+Ixx*Ixz*x(7,1)*x(8,1)-Ixz*Iyy*x(7,1)*x(8,1)+Ixz*Izz*x(7,1)*x(8,1)+Iyy*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz)-g*x(7,1)*cos(x(5,1))*sin(x(4,1)),g*(x(9,1)*cos(x(5,1))-x(7,1)*cos(x(4,1))*sin(x(5,1)))-((x(18,1)*cos(x(5,1))+x(16,1)*cos(x(6,1))*sin(x(5,1))+x(17,1)*sin(x(5,1))*sin(x(6,1)))*(Ixx*vtil(6,1)+Ixz*vtil(4,1)+Ixx^2*x(7,1)*x(8,1)+Ixz^2*x(7,1)*x(8,1)-Ixx*Ixz*x(8,1)*x(9,1)-Ixx*Iyy*x(7,1)*x(8,1)+Ixz*Iyy*x(8,1)*x(9,1)-Ixz*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz)-((Ixx*vtil(6,1)+Ixz*vtil(4,1))*(x(18,1)*cos(x(5,1))+x(16,1)*cos(x(6,1))*sin(x(5,1))+x(17,1)*sin(x(5,1))*sin(x(6,1))))/(Ixz^2-Ixx*Izz)-(cos(x(4,1))*(Ixz*vtil(6,1)+Izz*vtil(4,1))*(x(16,1)*cos(x(5,1))*cos(x(6,1))-x(18,1)*sin(x(5,1))+x(17,1)*cos(x(5,1))*sin(x(6,1))))/(Ixz^2-Ixx*Izz)-(cos(x(4,1))*(x(16,1)*cos(x(5,1))*cos(x(6,1))-x(18,1)*sin(x(5,1))+x(17,1)*cos(x(5,1))*sin(x(6,1)))*(Ixz*vtil(6,1)+Izz*vtil(4,1)-Ixz^2*x(8,1)*x(9,1)-Izz^2*x(8,1)*x(9,1)+Ixx*Ixz*x(7,1)*x(8,1)-Ixz*Iyy*x(7,1)*x(8,1)+Ixz*Izz*x(7,1)*x(8,1)+Iyy*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz),((x(17,1)*cos(x(5,1))*cos(x(6,1))-x(16,1)*cos(x(5,1))*sin(x(6,1)))*(Ixx*vtil(6,1)+Ixz*vtil(4,1)+Ixx^2*x(7,1)*x(8,1)+Ixz^2*x(7,1)*x(8,1)-Ixx*Ixz*x(8,1)*x(9,1)-Ixx*Iyy*x(7,1)*x(8,1)+Ixz*Iyy*x(8,1)*x(9,1)-Ixz*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz)-((x(16,1)*(cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(17,1)*(sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1))))*(Ixz*vtil(6,1)+Izz*vtil(4,1)-Ixz^2*x(8,1)*x(9,1)-Izz^2*x(8,1)*x(9,1)+Ixx*Ixz*x(7,1)*x(8,1)-Ixz*Iyy*x(7,1)*x(8,1)+Ixz*Izz*x(7,1)*x(8,1)+Iyy*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz)-((x(16,1)*(cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(17,1)*(sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1))))*(Ixz*vtil(6,1)+Izz*vtil(4,1)))/(Ixz^2-Ixx*Izz)+((x(17,1)*cos(x(5,1))*cos(x(6,1))-x(16,1)*cos(x(5,1))*sin(x(6,1)))*(Ixx*vtil(6,1)+Ixz*vtil(4,1)))/(Ixz^2-Ixx*Izz),g*cos(x(4,1))*cos(x(5,1))+(x(8,1)*(Ixx^2-Iyy*Ixx+Ixz^2)*(x(13,1)-x(18,1)*sin(x(5,1))+x(16,1)*cos(x(5,1))*cos(x(6,1))+x(17,1)*cos(x(5,1))*sin(x(6,1))))/(Ixz^2-Ixx*Izz)-(Ixz*x(8,1)*(Ixx-Iyy+Izz)*(x(15,1)+x(18,1)*cos(x(4,1))*cos(x(5,1))-x(17,1)*cos(x(6,1))*sin(x(4,1))+x(16,1)*sin(x(4,1))*sin(x(6,1))+x(16,1)*cos(x(4,1))*cos(x(6,1))*sin(x(5,1))+x(17,1)*cos(x(4,1))*sin(x(5,1))*sin(x(6,1))))/(Ixz^2-Ixx*Izz),((x(15,1)+x(16,1)*(sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1)))-x(17,1)*(cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(18,1)*cos(x(4,1))*cos(x(5,1)))*(Ixz^2*x(9,1)+Izz^2*x(9,1)-Ixx*Ixz*x(7,1)+Ixz*Iyy*x(7,1)-Ixz*Izz*x(7,1)-Iyy*Izz*x(9,1)))/(Ixz^2-Ixx*Izz)+((x(13,1)-x(18,1)*sin(x(5,1))+x(16,1)*cos(x(5,1))*cos(x(6,1))+x(17,1)*cos(x(5,1))*sin(x(6,1)))*(Ixx^2*x(7,1)+Ixz^2*x(7,1)-Ixx*Ixz*x(9,1)-Ixx*Iyy*x(7,1)+Ixz*Iyy*x(9,1)-Ixz*Izz*x(9,1)))/(Ixz^2-Ixx*Izz),g*sin(x(5,1))+(x(8,1)*(Ixz^2+Izz^2-Iyy*Izz)*(x(15,1)+x(18,1)*cos(x(4,1))*cos(x(5,1))-x(17,1)*cos(x(6,1))*sin(x(4,1))+x(16,1)*sin(x(4,1))*sin(x(6,1))+x(16,1)*cos(x(4,1))*cos(x(6,1))*sin(x(5,1))+x(17,1)*cos(x(4,1))*sin(x(5,1))*sin(x(6,1))))/(Ixz^2-Ixx*Izz)-(Ixz*x(8,1)*(Ixx-Iyy+Izz)*(x(13,1)-x(18,1)*sin(x(5,1))+x(16,1)*cos(x(5,1))*cos(x(6,1))+x(17,1)*cos(x(5,1))*sin(x(6,1))))/(Ixz^2-Ixx*Izz),0,0,0,(2*Ixx*vtil(6,1)+2*Ixz*vtil(4,1)+Ixx^2*x(7,1)*x(8,1)+Ixz^2*x(7,1)*x(8,1)-Ixx*Ixz*x(8,1)*x(9,1)-Ixx*Iyy*x(7,1)*x(8,1)+Ixz*Iyy*x(8,1)*x(9,1)-Ixz*Izz*x(8,1)*x(9,1))/(Ixz^2-Ixx*Izz),0,-(2*Ixz*vtil(6,1)+2*Izz*vtil(4,1)-Ixz^2*x(8,1)*x(9,1)-Izz^2*x(8,1)*x(9,1)+Ixx*Ixz*x(7,1)*x(8,1)-Ixz*Iyy*x(7,1)*x(8,1)+Ixz*Izz*x(7,1)*x(8,1)+Iyy*Izz*x(8,1)*x(9,1))/(Ixz^2-Ixx*Izz),(cos(x(5,1))*cos(x(6,1))*(Ixx*vtil(6,1)+Ixz*vtil(4,1)))/(Ixz^2-Ixx*Izz)-((sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1)))*(Ixz*vtil(6,1)+Izz*vtil(4,1)))/(Ixz^2-Ixx*Izz)-((sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1)))*(Ixz*vtil(6,1)+Izz*vtil(4,1)-Ixz^2*x(8,1)*x(9,1)-Izz^2*x(8,1)*x(9,1)+Ixx*Ixz*x(7,1)*x(8,1)-Ixz*Iyy*x(7,1)*x(8,1)+Ixz*Izz*x(7,1)*x(8,1)+Iyy*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz)+(cos(x(5,1))*cos(x(6,1))*(Ixx*vtil(6,1)+Ixz*vtil(4,1)+Ixx^2*x(7,1)*x(8,1)+Ixz^2*x(7,1)*x(8,1)-Ixx*Ixz*x(8,1)*x(9,1)-Ixx*Iyy*x(7,1)*x(8,1)+Ixz*Iyy*x(8,1)*x(9,1)-Ixz*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz),((cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1)))*(Ixz*vtil(6,1)+Izz*vtil(4,1)-Ixz^2*x(8,1)*x(9,1)-Izz^2*x(8,1)*x(9,1)+Ixx*Ixz*x(7,1)*x(8,1)-Ixz*Iyy*x(7,1)*x(8,1)+Ixz*Izz*x(7,1)*x(8,1)+Iyy*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz)+((cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1)))*(Ixz*vtil(6,1)+Izz*vtil(4,1)))/(Ixz^2-Ixx*Izz)+(cos(x(5,1))*sin(x(6,1))*(Ixx*vtil(6,1)+Ixz*vtil(4,1)))/(Ixz^2-Ixx*Izz)+(cos(x(5,1))*sin(x(6,1))*(Ixx*vtil(6,1)+Ixz*vtil(4,1)+Ixx^2*x(7,1)*x(8,1)+Ixz^2*x(7,1)*x(8,1)-Ixx*Ixz*x(8,1)*x(9,1)-Ixx*Iyy*x(7,1)*x(8,1)+Ixz*Iyy*x(8,1)*x(9,1)-Ixz*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz),-(sin(x(5,1))*(Ixx*vtil(6,1)+Ixz*vtil(4,1)))/(Ixz^2-Ixx*Izz)-(sin(x(5,1))*(Ixx*vtil(6,1)+Ixz*vtil(4,1)+Ixx^2*x(7,1)*x(8,1)+Ixz^2*x(7,1)*x(8,1)-Ixx*Ixz*x(8,1)*x(9,1)-Ixx*Iyy*x(7,1)*x(8,1)+Ixz*Iyy*x(8,1)*x(9,1)-Ixz*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz)-(cos(x(4,1))*cos(x(5,1))*(Ixz*vtil(6,1)+Izz*vtil(4,1)))/(Ixz^2-Ixx*Izz)-(cos(x(4,1))*cos(x(5,1))*(Ixz*vtil(6,1)+Izz*vtil(4,1)-Ixz^2*x(8,1)*x(9,1)-Izz^2*x(8,1)*x(9,1)+Ixx*Ixz*x(7,1)*x(8,1)-Ixz*Iyy*x(7,1)*x(8,1)+Ixz*Izz*x(7,1)*x(8,1)+Iyy*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz);
    0,0,0,((Ixz*vtil(6,1)+Izz*vtil(4,1))*(x(16,1)*(sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1)))-x(17,1)*(cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(18,1)*cos(x(4,1))*cos(x(5,1))))/(Ixz^2-Ixx*Izz)+((x(16,1)*(sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1)))-x(17,1)*(cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(18,1)*cos(x(4,1))*cos(x(5,1)))*(Ixz*vtil(6,1)+Izz*vtil(4,1)-Ixz^2*x(8,1)*x(9,1)-Izz^2*x(8,1)*x(9,1)+Ixx*Ixz*x(7,1)*x(8,1)-Ixz*Iyy*x(7,1)*x(8,1)+Ixz*Izz*x(7,1)*x(8,1)+Iyy*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz)-g*x(7,1)*cos(x(4,1))*cos(x(5,1)),(sin(x(4,1))*(x(16,1)*cos(x(5,1))*cos(x(6,1))-x(18,1)*sin(x(5,1))+x(17,1)*cos(x(5,1))*sin(x(6,1)))*(Ixz*vtil(6,1)+Izz*vtil(4,1)-Ixz^2*x(8,1)*x(9,1)-Izz^2*x(8,1)*x(9,1)+Ixx*Ixz*x(7,1)*x(8,1)-Ixz*Iyy*x(7,1)*x(8,1)+Ixz*Izz*x(7,1)*x(8,1)+Iyy*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz)-((x(18,1)*cos(x(5,1))+x(16,1)*cos(x(6,1))*sin(x(5,1))+x(17,1)*sin(x(5,1))*sin(x(6,1)))*(vtil(5,1)-Ixz*x(7,1)^2+Ixz*x(9,1)^2-Ixx*x(7,1)*x(9,1)+Izz*x(7,1)*x(9,1)))/Iyy-(vtil(5,1)*(x(18,1)*cos(x(5,1))+x(16,1)*cos(x(6,1))*sin(x(5,1))+x(17,1)*sin(x(5,1))*sin(x(6,1))))/Iyy-g*x(8,1)*cos(x(5,1))+g*x(7,1)*sin(x(4,1))*sin(x(5,1))+(sin(x(4,1))*(Ixz*vtil(6,1)+Izz*vtil(4,1))*(x(16,1)*cos(x(5,1))*cos(x(6,1))-x(18,1)*sin(x(5,1))+x(17,1)*cos(x(5,1))*sin(x(6,1))))/(Ixz^2-Ixx*Izz),(vtil(5,1)*(x(17,1)*cos(x(5,1))*cos(x(6,1))-x(16,1)*cos(x(5,1))*sin(x(6,1))))/Iyy-((x(16,1)*(cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(17,1)*(cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1))))*(Ixz*vtil(6,1)+Izz*vtil(4,1)-Ixz^2*x(8,1)*x(9,1)-Izz^2*x(8,1)*x(9,1)+Ixx*Ixz*x(7,1)*x(8,1)-Ixz*Iyy*x(7,1)*x(8,1)+Ixz*Izz*x(7,1)*x(8,1)+Iyy*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz)-((x(16,1)*(cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(17,1)*(cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1))))*(Ixz*vtil(6,1)+Izz*vtil(4,1)))/(Ixz^2-Ixx*Izz)+((x(17,1)*cos(x(5,1))*cos(x(6,1))-x(16,1)*cos(x(5,1))*sin(x(6,1)))*(vtil(5,1)-Ixz*x(7,1)^2+Ixz*x(9,1)^2-Ixx*x(7,1)*x(9,1)+Izz*x(7,1)*x(9,1)))/Iyy,(Ixz*x(8,1)*(Ixx-Iyy+Izz)*(x(14,1)+x(17,1)*cos(x(4,1))*cos(x(6,1))-x(16,1)*cos(x(4,1))*sin(x(6,1))+x(18,1)*cos(x(5,1))*sin(x(4,1))+x(16,1)*cos(x(6,1))*sin(x(4,1))*sin(x(5,1))+x(17,1)*sin(x(4,1))*sin(x(5,1))*sin(x(6,1))))/(Ixz^2-Ixx*Izz)-((Ixx*x(9,1)+2*Ixz*x(7,1)-Izz*x(9,1))*(x(13,1)-x(18,1)*sin(x(5,1))+x(16,1)*cos(x(5,1))*cos(x(6,1))+x(17,1)*cos(x(5,1))*sin(x(6,1))))/Iyy-g*cos(x(5,1))*sin(x(4,1)),-g*sin(x(5,1))-((x(14,1)-x(16,1)*(cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1)))+x(17,1)*(cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(18,1)*cos(x(5,1))*sin(x(4,1)))*(Ixz^2*x(9,1)+Izz^2*x(9,1)-Ixx*Ixz*x(7,1)+Ixz*Iyy*x(7,1)-Ixz*Izz*x(7,1)-Iyy*Izz*x(9,1)))/(Ixz^2-Ixx*Izz),((2*Ixz*x(9,1)-Ixx*x(7,1)+Izz*x(7,1))*(x(13,1)-x(18,1)*sin(x(5,1))+x(16,1)*cos(x(5,1))*cos(x(6,1))+x(17,1)*cos(x(5,1))*sin(x(6,1))))/Iyy-((x(8,1)*Ixz^2+x(8,1)*Izz^2-Iyy*x(8,1)*Izz)*(x(14,1)-x(16,1)*(cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1)))+x(17,1)*(cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(18,1)*cos(x(5,1))*sin(x(4,1))))/(Ixz^2-Ixx*Izz),0,0,0,(2*vtil(5,1)-Ixz*x(7,1)^2+Ixz*x(9,1)^2-Ixx*x(7,1)*x(9,1)+Izz*x(7,1)*x(9,1))/Iyy,(2*Ixz*vtil(6,1)+2*Izz*vtil(4,1)-Ixz^2*x(8,1)*x(9,1)-Izz^2*x(8,1)*x(9,1)+Ixx*Ixz*x(7,1)*x(8,1)-Ixz*Iyy*x(7,1)*x(8,1)+Ixz*Izz*x(7,1)*x(8,1)+Iyy*Izz*x(8,1)*x(9,1))/(Ixz^2-Ixx*Izz),0,(vtil(5,1)*cos(x(5,1))*cos(x(6,1)))/Iyy-((cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1)))*(Ixz*vtil(6,1)+Izz*vtil(4,1)))/(Ixz^2-Ixx*Izz)-((cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1)))*(Ixz*vtil(6,1)+Izz*vtil(4,1)-Ixz^2*x(8,1)*x(9,1)-Izz^2*x(8,1)*x(9,1)+Ixx*Ixz*x(7,1)*x(8,1)-Ixz*Iyy*x(7,1)*x(8,1)+Ixz*Izz*x(7,1)*x(8,1)+Iyy*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz)+(cos(x(5,1))*cos(x(6,1))*(vtil(5,1)-Ixz*x(7,1)^2+Ixz*x(9,1)^2-Ixx*x(7,1)*x(9,1)+Izz*x(7,1)*x(9,1)))/Iyy,((cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))*(Ixz*vtil(6,1)+Izz*vtil(4,1)-Ixz^2*x(8,1)*x(9,1)-Izz^2*x(8,1)*x(9,1)+Ixx*Ixz*x(7,1)*x(8,1)-Ixz*Iyy*x(7,1)*x(8,1)+Ixz*Izz*x(7,1)*x(8,1)+Iyy*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz)+((cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))*(Ixz*vtil(6,1)+Izz*vtil(4,1)))/(Ixz^2-Ixx*Izz)+(vtil(5,1)*cos(x(5,1))*sin(x(6,1)))/Iyy+(cos(x(5,1))*sin(x(6,1))*(vtil(5,1)-Ixz*x(7,1)^2+Ixz*x(9,1)^2-Ixx*x(7,1)*x(9,1)+Izz*x(7,1)*x(9,1)))/Iyy,(cos(x(5,1))*sin(x(4,1))*(Ixz*vtil(6,1)+Izz*vtil(4,1)))/(Ixz^2-Ixx*Izz)-(sin(x(5,1))*(vtil(5,1)-Ixz*x(7,1)^2+Ixz*x(9,1)^2-Ixx*x(7,1)*x(9,1)+Izz*x(7,1)*x(9,1)))/Iyy-(vtil(5,1)*sin(x(5,1)))/Iyy+(cos(x(5,1))*sin(x(4,1))*(Ixz*vtil(6,1)+Izz*vtil(4,1)-Ixz^2*x(8,1)*x(9,1)-Izz^2*x(8,1)*x(9,1)+Ixx*Ixz*x(7,1)*x(8,1)-Ixz*Iyy*x(7,1)*x(8,1)+Ixz*Izz*x(7,1)*x(8,1)+Iyy*Izz*x(8,1)*x(9,1)))/(Ixz^2-Ixx*Izz);
    0,0,0,0,vtil(9,1)*cos(x(5,1))+vtil(7,1)*cos(x(6,1))*sin(x(5,1))+vtil(8,1)*sin(x(5,1))*sin(x(6,1)),-cos(x(5,1))*(vtil(8,1)*cos(x(6,1))-vtil(7,1)*sin(x(6,1))),0,-x(15,1),x(14,1),1,0,0,0,x(9,1),-x(8,1),0,0,0;
    0,0,0,vtil(8,1)*(cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1)))-vtil(7,1)*(sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1)))-vtil(9,1)*cos(x(4,1))*cos(x(5,1)),-sin(x(4,1))*(vtil(7,1)*cos(x(5,1))*cos(x(6,1))-vtil(9,1)*sin(x(5,1))+vtil(8,1)*cos(x(5,1))*sin(x(6,1))),vtil(7,1)*(cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))+vtil(8,1)*(cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1))),x(15,1),0,-x(13,1),0,1,0,-x(9,1),0,x(7,1),0,0,0;
    0,0,0,vtil(8,1)*(cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))-vtil(7,1)*(cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1)))+vtil(9,1)*cos(x(5,1))*sin(x(4,1)),-cos(x(4,1))*(vtil(7,1)*cos(x(5,1))*cos(x(6,1))-vtil(9,1)*sin(x(5,1))+vtil(8,1)*cos(x(5,1))*sin(x(6,1))),-vtil(7,1)*(cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1)))-vtil(8,1)*(sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1))),-x(14,1),x(13,1),0,0,0,1,x(8,1),-x(7,1),0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
D = [0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0;
    0,0,0,0,0,0,0,0,0;
    0,0,0,-Izz/(Ixz^2-Ixx*Izz),0,-Ixz/(Ixz^2-Ixx*Izz),0,0,0;
    0,0,0,0,1/Iyy,0,0,0,0;
    0,0,0,-Ixz/(Ixz^2-Ixx*Izz),0,-Ixx/(Ixz^2-Ixx*Izz),0,0,0;
    1/m,0,0,-(2*Ixz*(x(14,1)-x(16,1)*(cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1)))+x(17,1)*(cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(18,1)*cos(x(5,1))*sin(x(4,1))))/(Ixz^2-Ixx*Izz),-(2*(x(15,1)+x(16,1)*(sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1)))-x(17,1)*(cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(18,1)*cos(x(4,1))*cos(x(5,1))))/Iyy,-(2*Ixx*(x(14,1)-x(16,1)*(cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1)))+x(17,1)*(cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(18,1)*cos(x(5,1))*sin(x(4,1))))/(Ixz^2-Ixx*Izz),0,0,0;
    0,1/m,0,(2*Ixz*(x(13,1)-x(18,1)*sin(x(5,1))+x(16,1)*cos(x(5,1))*cos(x(6,1))+x(17,1)*cos(x(5,1))*sin(x(6,1))))/(Ixz^2-Ixx*Izz)-(2*Izz*(x(15,1)+x(16,1)*(sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1)))-x(17,1)*(cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(18,1)*cos(x(4,1))*cos(x(5,1))))/(Ixz^2-Ixx*Izz),0,(2*Ixx*(x(13,1)-x(18,1)*sin(x(5,1))+x(16,1)*cos(x(5,1))*cos(x(6,1))+x(17,1)*cos(x(5,1))*sin(x(6,1))))/(Ixz^2-Ixx*Izz)-(2*Ixz*(x(15,1)+x(16,1)*(sin(x(4,1))*sin(x(6,1))+cos(x(4,1))*cos(x(6,1))*sin(x(5,1)))-x(17,1)*(cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(18,1)*cos(x(4,1))*cos(x(5,1))))/(Ixz^2-Ixx*Izz),0,0,0;
    0,0,1/m,(2*Izz*(x(14,1)-x(16,1)*(cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1)))+x(17,1)*(cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(18,1)*cos(x(5,1))*sin(x(4,1))))/(Ixz^2-Ixx*Izz),(2*(x(13,1)-x(18,1)*sin(x(5,1))+x(16,1)*cos(x(5,1))*cos(x(6,1))+x(17,1)*cos(x(5,1))*sin(x(6,1))))/Iyy,(2*Ixz*(x(14,1)-x(16,1)*(cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1)))+x(17,1)*(cos(x(4,1))*cos(x(6,1))+sin(x(4,1))*sin(x(5,1))*sin(x(6,1)))+x(18,1)*cos(x(5,1))*sin(x(4,1))))/(Ixz^2-Ixx*Izz),0,0,0;
    0,0,0,0,0,0,-cos(x(5,1))*cos(x(6,1)),-cos(x(5,1))*sin(x(6,1)),sin(x(5,1));
    0,0,0,0,0,0,cos(x(4,1))*sin(x(6,1))-cos(x(6,1))*sin(x(4,1))*sin(x(5,1)),-cos(x(4,1))*cos(x(6,1))-sin(x(4,1))*sin(x(5,1))*sin(x(6,1)),-cos(x(5,1))*sin(x(4,1));
    0,0,0,0,0,0,-sin(x(4,1))*sin(x(6,1))-cos(x(4,1))*cos(x(6,1))*sin(x(5,1)),cos(x(6,1))*sin(x(4,1))-cos(x(4,1))*sin(x(5,1))*sin(x(6,1)),-cos(x(4,1))*cos(x(5,1));
    0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,1];

end