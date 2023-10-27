function [y,H] = measmodel_windplusjerk(t,x,u,dervflag,params)
%
%  Copyright (c) 2023 Jeremy W. Hopwood.  All rights reserved.
% 
%  This function is a template for specifying the measurement model,
%
%       Continuous-time:    y(t) = h(t,x,u)                         (1)
%         Discrete-time:    y(k) = h(k,x(k),u(k))                   (2)
%    Continous-discrete:   y(tk) = h(tk,x(tk),u(tk))                (3)
%
%  It is assumed for most filtering problems that measurement noise, w, is
%  additive such that
%
%                           z(k) = h(k,x(k),u(k)) + w(k)            (4)
%
%  and thus no information about w is required in this function. 
%  h(t,x,u), and its derivative with respect to x, H = dh/dx. This function
%  is for use in either static estimation problems such as nonlinear
%  least-squares or dyamic estimation problems such as extended Kalman
%  filtering.
%
%  Inputs:
%
%   t           The time at which h is evaluated. If the system is
%               discrete, t=k is the sample number, where the initial
%               condition occurs at k=0.
%
%   x           The nx x 1 state vector at time t (or sample k).
%
%   u           The nu x 1 control vector at time t (or sample k).
%
%   dervflag    A flag that determines whether (dervflag = 1) or not
%               (dervflag = 0) the partial derivative dh/dx must be
%               calculated. If dervflag = 0, then this output will be
%               an empty array.
%
%   params      Any data type that is used to pass through constant
%               parameters to the dynamics model.
%  
%  Outputs:
%
%   y           The nz x 1 output vector at time t or sample k.
%
%   H           The nz x nx Jacobian matrix of h with respect to x.
%               This output is needed  when performaing extended Kalman
%               filtering or Gauss-Newton estimaton, for example.
%

% Determine the dimensions of x and y and compute the outputs.
nx = 18;
ny = 12;

% Get necessary constants
% g = params.g;

% Compute the rotation matrix
% phi = x(4,1);
% theta = x(5,1);
% psi = x(6,1);
% R1 = [1,0,0;0,cos(phi),-sin(phi);0,sin(phi),cos(phi)];
% R2 = [cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)];
% R3 = [cos(psi),-sin(psi),0;sin(psi),cos(psi),0;0,0,1];
% R_IB = R3*R2*R1;
% e3 = [0;0;1];

% Assemble output
y = x(1:ny,1);
% y(10:12,1) = y(10:12,1) + g*R_IB'*e3; % add gravity to accelerometer

% Return if neither first derivatives nor second derivatives need to be
% calculated.
if dervflag == 0
    H = [];
    return
end

% Calculate the first derivative.
H = [eye(ny), zeros(ny,nx-ny)];
% H(10,5) = -g*cos(theta);
% H(11,4) = g*cos(phi)*cos(theta);
% H(11,5) = -g*sin(phi)*sin(theta);
% H(12,4) = -g*cos(theta)*sin(phi);
% H(12,5) = -g*cos(phi)*sin(theta);

end