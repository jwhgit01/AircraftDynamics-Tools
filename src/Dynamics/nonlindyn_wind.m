function [xdot,A,D] = nonlindyn_wind(t,x,u,vtil,dervflag,params)
%nonlindyn
%
% Copyright (c) 2023 Jeremy W. Hopwood. All rights reserved.
%
% This is a function defines the aircraft dynamics in a frozen,
% uniform wind field in the form required by the Estimation-Tools
% toolbox. It is a model for thecontinuous-time, nonlinear model of
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
%               b                  The drift vector field for the wind
%                                  model. A function handle that is a
%                                  function of x.
%
%               D                  The diffusion matrix function for the
%                                  wind model. A function handle that is a
%                                  function of x.
%
%               Sigma              A 3x3 constant matrix. The infinitesimal
%                                  covariance of the driving Weiner process
%                                  is Q = Sigma'*Sigma. i.e. Sigma=chol(Q)
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
rho = params.rho;
n_us = params.AerodynamicModel.NumberOfStates;

% Parse the state vector and compute attitude kinematics
x_us = zeros(0,1);
if strcmp(params.AttitudeParameterization,'EulerAngles')
    Theta = x(4:6,1);
    omega = x(7:9,1);
    vbr = x(10:12,1);
    w = x(13:15,1);
    if n_us > 0
        x_us = x(16:15+n_us,1);
    end
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
elseif strcmp(params.AttitudeParameterization,'Quaternion')
    quat = x(4:7,1);
    omega = x(8:10,1);
    vbr = x(11:13,1);
    w = x(14:16,1);
    if n_us > 0
        x_us = x(17:16+n_us,1);
    end
    R_IB = quat2rotmat(quat);
    Omega = [-cpem(omega),omega;-omega.',0];
    attitudeDyn = 0.5*Omega*quat;
elseif strcmp(params.AttitudeParameterization,'CompassTilt')
    lambda = x(4:6,1);
    zeta = x(7:9,1);
    lambda = lambda/norm(lambda);
    zeta = zeta/norm(zeta);
    omega = x(10:12,1);
    vbr = x(13:15,1);
    w = x(16:18,1);
    if n_us > 0
        x_us = x(19:18+n_us,1);
    end
    R_IB = [lambda, cpem(zeta)*lambda, zeta].';
    attitudeDyn = [cross(lambda,omega); cross(zeta,omega)];
else
    error('Unknown attitude parameterization');
end

% translational kinematics
NEDdot = R_IB*vbr + w;

% Aerodynamic forces and moments
Parameters = params.ParameterEstimates;
Constants.rho = rho;
Constants.g = g;
[F,M] = params.AerodynamicModel.Model(x_us,vbr,omega,u,Parameters,Constants);
            
% Dynamic equations of motion
e3 = [0;0;1];
vbrdot = cross(vbr,omega) + g*R_IB'*e3 + F/m;
omegadot = I\(cross(I*omega,omega) + M);

% Wind dynamics
if isempty(vtil)
    vtil = zeros(3,1);
end
% wdot = params.b(x) + params.D(x)*params.Sigma*vtil;
wdot = params.b(x) + params.D(x)*vtil;

% Rigid body in wind state derivative
x_rbdot = [NEDdot;attitudeDyn;omegadot;vbrdot;wdot];

% Unsteady aerodynamics (TODO)
% x_usdot = params.AerodynamicModel.AugmentedDynamics(x_us,vbr,omega,u,Parameters,Constants);

% Assemble state derivative
xdot = x_rbdot; % [x_rbdot; x_usdot] (TODO)

% If the Jacobians are not needed, set to empty arrays and return
if dervflag == 0
    A = [];
    D = [];
    return
end

% Compute the Jacobains of the aerodynamic model
[dfdvr,dfdomega,~,~,~] = params.AerodynamicModel.Jacobian([],vbr,omega,u,Parameters,Constants);

% Compute the Jacobian of the rigid body dynamics (including wind states)
if strcmp(params.AttitudeParameterization,'EulerAngles')
    error('TODO') % A_rb = 
elseif strcmp(params.AttitudeParameterization,'Quaternion')
    error('TODO') % A_rb = 
elseif strcmp(params.AttitudeParameterization,'CompassTilt')
    Ixx = I(1,1);
    Iyy = I(2,2);
    Izz = I(3,3);
    Ixz = -I(1,3);
    lambda1 = x(4,1);
    lambda2 = x(5,1);
    lambda3 = x(6,1);
    zeta1 = x(7,1);
    zeta2 = x(8,1);
    zeta3 = x(9,1);
    p = x(10,1);
    q = x(11,1);
    r = x(12,1);
    ur = x(13,1);
    vr = x(14,1);
    wr = x(15,1);
    A_rb = [0,0,0,ur,vr,wr,0,0,0,0,0,0,lambda1,lambda2,lambda3,1,0,0;
            0,0,0,vr*zeta3-wr*zeta2,wr*zeta1-ur*zeta3,ur*zeta2-vr*zeta1,lambda2*wr-lambda3*vr,lambda3*ur-lambda1*wr,lambda1*vr-lambda2*ur,0,0,0,lambda3*zeta2-lambda2*zeta3,lambda1*zeta3-lambda3*zeta1,lambda2*zeta1-lambda1*zeta2,0,1,0;
            0,0,0,0,0,0,ur,vr,wr,0,0,0,zeta1,zeta2,zeta3,0,0,1;
            0,0,0,0,r,-q,0,0,0,0,-lambda3,lambda2,0,0,0,0,0,0;
            0,0,0,-r,0,p,0,0,0,lambda3,0,-lambda1,0,0,0,0,0,0;
            0,0,0,q,-p,0,0,0,0,-lambda2,lambda1,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,r,-q,0,-zeta3,zeta2,0,0,0,0,0,0;
            0,0,0,0,0,0,-r,0,p,zeta3,0,-zeta1,0,0,0,0,0,0;
            0,0,0,0,0,0,q,-p,0,-zeta2,zeta1,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0,0,-(q*(Ixx*Ixz-Ixz*Iyy+Ixz*Izz))/(Ixz^2-Ixx*Izz),(Ixz^2*r+Izz^2*r-Ixx*Ixz*p+Ixz*Iyy*p-Ixz*Izz*p-Iyy*Izz*r)/(Ixz^2-Ixx*Izz),(q*(Ixz^2+Izz^2-Iyy*Izz))/(Ixz^2-Ixx*Izz),0,0,0,0,0,0;
            0,0,0,0,0,0,0,0,0,-(2*Ixz*p+Ixx*r-Izz*r)/Iyy,0,(Izz*p-Ixx*p+2*Ixz*r)/Iyy,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0,0,-(q*(Ixx^2-Iyy*Ixx+Ixz^2))/(Ixz^2-Ixx*Izz),-(Ixx^2*p+Ixz^2*p-Ixx*Iyy*p-Ixx*Ixz*r+Ixz*Iyy*r-Ixz*Izz*r)/(Ixz^2-Ixx*Izz),(q*(Ixx*Ixz-Ixz*Iyy+Ixz*Izz))/(Ixz^2-Ixx*Izz),0,0,0,0,0,0;
            0,0,0,0,0,0,g,0,0,0,-wr,vr,0,r,-q,0,0,0;
            0,0,0,0,0,0,0,g,0,wr,0,-ur,-r,0,p,0,0,0;
            0,0,0,0,0,0,0,0,g,-vr,ur,0,q,-p,0,0,0,0;
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
            0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

else
    error('Unknown attitude parameterization');
end

% Construct the Jacobians
nv = 3;
nx = size(xdot,1);
n_att = size(attitudeDyn,1);
A_aero = [zeros(3+n_att,nx); 
          zeros(6,3+n_att),dfdomega,dfdvr,zeros(6,3);
          zeros(3,nx)];
A = A_rb + A_aero;
D = [zeros(3+n_att,nv); dfdvr; eye(3)];

end