classdef DrydenTurbulence < TurbulenceModel
%DrydenTurbulence 
%
% Copyright (c) 2023 Jeremy W. Hopwood. All rights reserved.
%
% This class ...
%

properties
    Altitude
    Specification
    Windspeed
    Airspeed
    Span
    StateSpaceModel
    WhiteNoiseFunction
    Simulation
end
    
methods

    function obj = DrydenTurbulence(Altitude,Specification,Windspeed,Airspeed,Span,WhiteNoiseFunction)
        %DrydenTurbulence
        obj.Altitude = Altitude;
        obj.Specification = Specification;
        obj.Windspeed = Windspeed;
        obj.Airspeed = Airspeed;
        obj.Span = Span;
        obj.WhiteNoiseFunction = WhiteNoiseFunction;

        % Create transfer function and realize it.
        H_tf = obj.DrydenTransferFunction;
        obj.StateSpaceModel = ss(H_tf);
        obj.NumberOfStates = size(obj.StateSpaceModel.A,1);
    end

    function [dw,dgradW] = Output(obj,t,x_env,x_rb)

        % Unit variance white Gaussian noise
        nu = obj.WhiteNoiseFunction(t);

        % Shaping filter output (body frame wind)
        y = obj.StateSpaceModel.C*x_env + obj.StateSpaceModel.D*nu;     
        dwb = y(1:3,1);
        domegawb = y(4:6,1);

        % Compute the attitude rotation matrix
        if strcmp(obj.Simulation.AttitudeParameterization,'EulerAngles')
            Theta = x_rb(4:6,1);
            phi = Theta(1);
            theta = Theta(2);
            psi = Theta(3);
            R1 = [1,0,0;0,cos(phi),-sin(phi);0,sin(phi),cos(phi)];
            R2 = [cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)];
            R3 = [cos(psi),-sin(psi),0;sin(psi),cos(psi),0;0,0,1];
            R_IB = R3*R2*R1;
        elseif strcmp(obj.Simulation.AttitudeParameterization,'Quaternion')
            quat = x_rb(4:7,1);
            R_IB = quat2rotmat(quat);
        elseif strcmp(obj.Simulation.AttitudeParameterization,'CompassTilt')
            lambda = x(4:6,1);
            zeta = x_rb(7:9,1);
            R_IB = [lambda, cpem(zeta)*lambda, zeta].';
        else
            error('Unknown attitude parameterization');
        end

        % Transform the wind and its grdient to the inertial frame.
        dw = R_IB*dwb; 
        Phi = zeros(3,3);
        Phi(3,2) = domegawb(1);
        Phi(3,1) = -domegawb(2);
        Phi(2,1) = domegawb(3);
        dgradW = R_IB*Phi*R_IB';

    end
    
    function x_envdot = Dynamics(obj,t,x_env,~)

        % Unit variance white Gaussian noise
        nu = obj.WhiteNoiseFunction(t);

        % Shaping filter dynamics
        x_envdot = obj.StateSpaceModel.A*x_env + obj.StateSpaceModel.B*nu;

    end

    function H = DrydenTransferFunction(obj)
    %DrydenTransferFunctions 
    % Source: https://www.mathworks.com/help/aeroblks/drydenwindturbulencemodelcontinuous.html
    
    % Determine characteristic lengths, L for low altitude turbulence (<1000ft)
    if strcmp(obj.Specification,'MIL-F-8785C')
        Lu = obj.Altitude/((0.177+0.000823*obj.Altitude)^1.2);
        Lv = Lu;
        Lw = obj.Altitude;
    elseif strcmp(obj.Specification,'MIL-HDBK-1797B')
        Lu = obj.Altitude/((0.177+0.000823*obj.Altitude)^1.2);
        Lv = Lu/2;
        Lw = obj.Altitude/2;
    else
        error('Invalid Specification');
    end
    
    % Turbulence standard deviations using low altitude model
    sigma_w = 0.1*obj.Windspeed;
    sigma_u = sigma_w/((0.177+0.000823*obj.Altitude)^1.2);
    sigma_v = sigma_u;
    
    % continuous time transfer functions
    if strcmp(obj.Specification,'MIL-F-8785C')
        
        % Longitudinal
        H(1,1) = tf(sigma_u*sqrt(2*Lu/(pi*obj.Airspeed)),[Lu/obj.Airspeed,1]);
        H(4,1) = tf(sigma_w*sqrt(0.8/obj.Airspeed)*(pi/(4*obj.Span))^(1/6),Lw^(1/3)*[4*obj.Span/(pi*obj.Airspeed),1]);
    
        % Lateral
        H(2,1) = tf(sigma_v*sqrt(Lv/(pi*obj.Airspeed))*[sqrt(3)*Lv/obj.Airspeed,1],[(Lv/obj.Airspeed)^2,2*Lv/obj.Airspeed,1]);
        H(5,1) = tf([1/obj.Airspeed,0],[3*obj.Span/(pi*obj.Airspeed),1])*H(2,1);
    
        % Vertical
        H(3,1) = tf(sigma_w*sqrt(Lw/(pi*obj.Airspeed))*[sqrt(3)*Lw/obj.Airspeed,1],[(Lw/obj.Airspeed)^2,2*Lw/obj.Airspeed,1]);
        H(6,1) = tf([1/obj.Airspeed,0],[4*obj.Span/(pi*obj.Airspeed),1])*H(3,1);
    
    elseif strcmp(obj.Specification,'MIL-HDBK-1797B')
        
        % Longitudinal
        H(1,1) = tf(sigma_u*sqrt(2*Lu/(pi*obj.Airspeed)),[Lu/obj.Airspeed,1]); % same as MIL-F-8785C
        H(4,1) = tf(sigma_w*sqrt(0.8/obj.Airspeed)*(pi/(4*obj.Span))^(1/6),(2*Lw)^(1/3)*[4*obj.Span/(pi*obj.Airspeed),1]);
    
        % Lateral
        H(2,1) = tf(sigma_v*sqrt(2*Lv/(pi*obj.Airspeed))*[2*sqrt(3)*Lv/obj.Airspeed,1],[(2*Lv/obj.Airspeed)^2,4*Lv/obj.Airspeed,1]);
        H(5,1) = tf([1/obj.Airspeed,0],[3*obj.Span/(pi*obj.Airspeed),1])*H(2,1); % same form as MIL-F-8785C
    
        % Vertical
        H(3,1) = tf(sigma_w*sqrt(2*Lw/(pi*obj.Airspeed))*[2*sqrt(3)*Lw/obj.Airspeed,1],[(2*Lw/obj.Airspeed)^2,4*Lw/obj.Airspeed,1]);
        H(6,1) = tf([1/obj.Airspeed,0],[4*obj.Span/(pi*obj.Airspeed),1])*H(3,1); % same form as MIL-F-8785C
    
    else
        error('Invalid Specification');
    end
    
    end

end % public Methods

end