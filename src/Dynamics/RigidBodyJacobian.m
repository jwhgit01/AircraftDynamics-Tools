function A = RigidBodyJacobian(x,I,g,AttitudeParam)
%RigidBodyJacobian
%
% Copyright (c) 2023 Jeremy W. Hopwood. All rights reserved.
%
% This function ...
%
% Inputs:
%
%   in1     The first input...
%  
% Outputs:
%
%   out1    The first output...
%

% Parse state vector and compute attitude kinematics
if strcmp(AttitudeParam,'EulerAngles')
    
    % aircraft states
    phi = x(4,1);
    theta = x(5,1);
    psi = x(6,1);
    u = x(7,1);
    v = x(8,1);
    w = x(9,1);
    p = x(10,1);
    q = x(11,1);
    r = x(12,1);
    
    % Jacobian matrix
    A =[0,0,0,v*(sin(phi)*sin(psi)+cos(phi)*cos(psi)*sin(theta))+w*(cos(phi)*sin(psi)-cos(psi)*sin(phi)*sin(theta)),cos(psi)*(w*cos(phi)*cos(theta)-u*sin(theta)+v*cos(theta)*sin(phi)),w*(cos(psi)*sin(phi)-cos(phi)*sin(psi)*sin(theta))-v*(cos(phi)*cos(psi)+sin(phi)*sin(psi)*sin(theta))-u*cos(theta)*sin(psi),cos(psi)*cos(theta),cos(psi)*sin(phi)*sin(theta)-cos(phi)*sin(psi),sin(phi)*sin(psi)+cos(phi)*cos(psi)*sin(theta),0,0,0;
        0,0,0,-v*(cos(psi)*sin(phi)-cos(phi)*sin(psi)*sin(theta))-w*(cos(phi)*cos(psi)+sin(phi)*sin(psi)*sin(theta)),sin(psi)*(w*cos(phi)*cos(theta)-u*sin(theta)+v*cos(theta)*sin(phi)),w*(sin(phi)*sin(psi)+cos(phi)*cos(psi)*sin(theta))-v*(cos(phi)*sin(psi)-cos(psi)*sin(phi)*sin(theta))+u*cos(psi)*cos(theta),cos(theta)*sin(psi),cos(phi)*cos(psi)+sin(phi)*sin(psi)*sin(theta),cos(phi)*sin(psi)*sin(theta)-cos(psi)*sin(phi),0,0,0;
        0,0,0,cos(theta)*(v*cos(phi)-w*sin(phi)),-u*cos(theta)-w*cos(phi)*sin(theta)-v*sin(phi)*sin(theta),0,-sin(theta),cos(theta)*sin(phi),cos(phi)*cos(theta),0,0,0;
        0,0,0,(sin(theta)*(q*cos(phi)-r*sin(phi)))/cos(theta),(r*cos(phi)+q*sin(phi))/cos(theta)^2,0,0,0,0,1,sin(phi)*tan(theta),cos(phi)*tan(theta);
        0,0,0,-r*cos(phi)-q*sin(phi),0,0,0,0,0,0,cos(phi),-sin(phi);
        0,0,0,(q*cos(phi)-r*sin(phi))/cos(theta),(sin(theta)*(r*(2*sin(phi/2)^2-1)-q*sin(phi)))/(sin(theta)^2-1),0,0,0,0,0,sin(phi)/cos(theta),cos(phi)/cos(theta);
        0,0,0,0,-g*cos(theta),0,0,r,-q,0,-w,v;
        0,0,0,g*cos(phi)*cos(theta),-g*sin(phi)*sin(theta),0,-r,0,p,w,0,-u;
        0,0,0,-g*cos(theta)*sin(phi),-g*cos(phi)*sin(theta),0,q,-p,0,-v,u,0;
        0,0,0,0,0,0,0,0,0,-((-I(1,3))*q*(I(1,1)-I(2,2)+I(3,3)))/((-I(1,3))^2-I(1,1)*I(3,3)),((-I(1,3))^2*r+I(3,3)^2*r-I(1,1)*(-I(1,3))*p+(-I(1,3))*I(2,2)*p-(-I(1,3))*I(3,3)*p-I(2,2)*I(3,3)*r)/((-I(1,3))^2-I(1,1)*I(3,3)),(q*((-I(1,3))^2+I(3,3)^2-I(2,2)*I(3,3)))/((-I(1,3))^2-I(1,1)*I(3,3));
        0,0,0,0,0,0,0,0,0,-(2*(-I(1,3))*p+I(1,1)*r-I(3,3)*r)/I(2,2),0,(I(3,3)*p-I(1,1)*p+2*(-I(1,3))*r)/I(2,2);
        0,0,0,0,0,0,0,0,0,-(q*(I(1,1)^2-I(2,2)*I(1,1)+(-I(1,3))^2))/((-I(1,3))^2-I(1,1)*I(3,3)),-(I(1,1)^2*p+(-I(1,3))^2*p-I(1,1)*I(2,2)*p-I(1,1)*(-I(1,3))*r+(-I(1,3))*I(2,2)*r-(-I(1,3))*I(3,3)*r)/((-I(1,3))^2-I(1,1)*I(3,3)),((-I(1,3))*q*(I(1,1)-I(2,2)+I(3,3)))/((-I(1,3))^2-I(1,1)*I(3,3))];

elseif strcmp(AttitudeParam,'Quaternion')
    error('TODO')
elseif strcmp(AttitudeParam,'CompassTilt')
    error('TODO')
else
    error('Unknown attitude parameterization');
end

end