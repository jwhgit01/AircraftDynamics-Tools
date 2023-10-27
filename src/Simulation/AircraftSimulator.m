classdef AircraftSimulator < handle
%AircraftSimulator 
%
% Copyright (c) 2023 Jeremy W. Hopwood. All rights reserved.
%
% This class ...
%


properties
    Name
    TimeStep
    FinalTime
    InitialCondition
    Aircraft
    AerodynamicModel
    AerodynamicParameters % ParameterEstimates Object or array
    ActuatorModel
    EnviromentModel
    SensorModel % 'IMU+GPS' or 'INS'
    AttitudeParameterization % 'EulerAngles','Quaternion', or 'CompassTilt'
    INS % insSensor object (Navigation Toolbox)
    GNSS % gnssSensor object (Navigation Toolbox)
    IMU % imuSensor object (Navigation Toolbox)
    ControlLaw % ControlLaw object
    Observer % Observer object
    ReferenceLocation % Origin of NED frame in the LLA frame [lat;lon;alt]
    Dimensions
end

properties (Access = protected)
    Acceleration = zeros(3,1);
    SimulationInitialized = false;
end

methods

    function obj = AircraftSimulator(Name)
        %AircraftSimulator Construct an instance of this class with
        % default properties.

        % Unique ID
        obj.Name = Name;

        % Default time vector
        obj.FinalTime = 30;
        obj.TimeStep = 0.01;

        % Default NED origin position
        obj.ReferenceLocation = [37.19711;-80.57810;528]; % KEAS

        % Create the default enviroment
        obj.EnviromentModel = EnviromentModel();

    end

    function tspan = Time(obj)
        dt = obj.TimeStep;
        T = obj.FinalTime;
        tspan = (0:dt:T).';
    end

    function Initialize(obj)

        % Numbers of rigid body dynamics states
        if strcmp(obj.AttitudeParameterization,'EulerAngles')
            n_rb = 12;
        elseif strcmp(obj.AttitudeParameterization,'Quaternion')
            n_rb = 13;
        elseif strcmp(obj.AttitudeParameterization,'CompassTilt')
            n_rb = 15;
        else
            error('Unknown attitude parameterization');
        end

        % Number of unsteady/quasi-steady states
        n_us = obj.AerodynamicModel.NumberOfStates;

        % Number of actuator states
        if isempty(obj.ActuatorModel)
            n_act = 0;
        else
            n_act = 0;
            for ii = 1:length(obj.ActuatorModel)
                n_act = n_act + obj.ActuatorModel(ii).NumberOfStates;
            end
        end
        
        % Number of enviroment states
        if isempty(obj.EnviromentModel.TurbulenceModel)
            n_env = 0;
        else
            n_env = obj.EnviromentModel.TurbulenceModel.NumberOfStates;
        end

        % Number of observer states
        if isempty(obj.Observer)
            n_obsv = 0;
        else
            n_obsv = obj.Observer.NumberOfStates;
        end

        % Number of controller states
        if isempty(obj.ControlLaw)
            n_ctrl = 0;
        else
            n_ctrl = obj.ControlLaw.NumberOfStates;
        end

        % Save dimensions in Dimensions property
        obj.Dimensions.n_rb = n_rb;
        obj.Dimensions.n_us = n_us;
        obj.Dimensions.n_act = n_act;
        obj.Dimensions.n_env = n_env;
        obj.Dimensions.n_obsv = n_obsv;
        obj.Dimensions.n_ctrl = n_ctrl;

        % Successful
        obj.SimulationInitialized = true;

    end

    function [x_rb,x_us,x_act,x_env,x_obsv,x_ctrl] = State(obj,x)

        % Get dimensions
        n_rb = obj.Dimensions.n_rb;
        n_us = obj.Dimensions.n_us;
        n_act = obj.Dimensions.n_act;
        n_env = obj.Dimensions.n_env;
        n_obsv = obj.Dimensions.n_obsv;
        n_ctrl = obj.Dimensions.n_ctrl;

        % If x is given as a vector,
        if size(x,2) == 1
            x_rb = x(1:n_rb,:);
            x_us = x(n_rb+1:n_rb+n_us,:);
            x_act = x(n_rb+n_us+1:n_rb+n_us+n_act,:);
            x_env = x(n_rb+n_us+n_act+1:n_rb+n_us+n_act+n_env,:);
            x_obsv = x(n_rb+n_us+n_act+n_env+1:n_rb+n_us+n_act+n_env+n_obsv,:);
            x_ctrl = x(n_rb+n_us+n_act+n_env+n_obsv+1:n_rb+n_us+n_act+n_env+n_obsv+n_ctrl,:);
            return
        end

        % If x is given as a time history,
        x_rb = x(:,1:n_rb);
        x_us = x(:,n_rb+1:n_rb+n_us);
        x_act = x(:,n_rb+n_us+1:n_rb+n_us+n_act);
        x_env = x(:,n_rb+n_us+n_act+1:n_rb+n_us+n_act+n_env);
        x_obsv = x(:,n_rb+n_us+n_act+n_env+1:n_rb+n_us+n_act+n_env+n_obsv);
        x_ctrl = x(:,n_rb+n_us+n_act+n_env+n_obsv+1:n_rb+n_us+n_act+n_env+n_obsv+n_ctrl);


    end

    function xdot = Dynamics(obj,t,x)

        % Check if the simulation was initialized
        if ~obj.SimulationInitialized
            error('Simulation not initialized. Run sim.Initialize')
        end

        % Parse state vector
        [x_rb,x_us,x_act,x_env,x_obsv,x_ctrl] = obj.State(x);

        % Position in local NED frame
        % NED = x_rb(1:3,1);

        % Parse rigid body state vector and compute attitude kinematics
        if strcmp(obj.AttitudeParameterization,'EulerAngles')
            Theta = x_rb(4:6,1);
            vb = x_rb(7:9,1);
            omega = x_rb(10:12,1);
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
        elseif strcmp(obj.AttitudeParameterization,'Quaternion')
            quat = x_rb(4:7,1);
            vb = x_rb(8:10,1);
            omega = x_rb(11:13,1);
            R_IB = quat2rotmat(quat);
            Omega = [-cpem(omega),omega;-omega.',0];
            attitudeDyn = 0.5*Omega*quat;
        elseif strcmp(obj.AttitudeParameterization,'CompassTilt')
            lambda = x_rb(4:6,1);
            zeta = x_rb(7:9,1);
            vb = x_rb(10:12,1);
            omega = x_rb(13:15,1);
            R_IB = [lambda, cpem(zeta)*lambda, zeta].';
            attitudeDyn = [cross(lambda,omega); cross(zeta,omega)];
        else
            error('Unknown attitude parameterization');
        end
         
        % translational kinematics
        NEDdot = R_IB*vb;

        % Wind velocity and gradient at current time and position
        [w,gradW] = obj.EnviromentModel.Output(t,x_env,x_rb);

        % Relative body velocity and angular velocity
        vb_r = vb - R_IB'*w;
        Phi = R_IB'*gradW*R_IB;
        omega_w = [Phi(3,2);-Phi(3,1);Phi(2,1)];
        omega_r = omega - omega_w;

        % Sensor outputs
%         if strcmp(obj.SensorModel,'INS')
%             motion.Orientation = quaternion(quat.');
%             motion.Position = NED.';
%             motion.Velocity = vb.';
%             motion.Acceleration = obj.Acceleration.'; % non-causal
%             motion.AngularVelocity = omega.';
%             y = obj.INS(motion);
%         elseif strcmp(obj.SensorModel,'IMU+GPS')
%             error('TODO')
%         else
%             error('Invalid SensorModel')
%         end
        y = [];

        % Controller outputs
        if isempty(obj.ControlLaw)
            u = [];
        else
            u = obj.ControlLaw.Output(t,x_ctrl,x_rb,x_obsv,y);
        end

        % Actuator outputs
        if isempty(obj.ActuatorModel)
            delta = u;
        else
            delta = obj.ActuatorModel.Output(x_act,u);
        end

        % Envirmental Constants
        rho = obj.EnviromentModel.AirDensity;
        g = obj.EnviromentModel.Gravity;

        % Aircraft Constants
        m = obj.Aircraft.Mass;
        I = obj.Aircraft.MomentOfInertia;
       
        % Aerodynamic forces and moments
        Parameters = obj.AerodynamicParameters;
        Constants.rho = rho;
        Constants.g = g;
        [F,M] = obj.AerodynamicModel.Model(x_us,vb_r,omega_r,delta,Parameters,Constants);
                    
        % Dynamic equations of motion
        e3 = [0;0;1];
        vbdot = cross(vb,omega) + g*R_IB'*e3 + F/m;
        omegadot = I\(cross(I*omega,omega) + M);

        % Store acceleration
%         obj.Acceleration = vbdot + cross(omega,vb);
                    
        % Rigid body state derivative
        x_rbdot = [NEDdot;attitudeDyn;vbdot;omegadot];

        % Unsteady aerodynamics
        x_usdot = obj.AerodynamicModel.AugmentedDynamics(x_us,vb_r,omega_r,delta,Parameters,Constants);

        % Actuator dynamics
        x_actdot = zeros(0,1);
        if ~isempty(obj.ActuatorModel)
            for ii = 1:length(obj.ActuatorModel)
                x_actdot = [x_actdot;obj.ActuatorModel(ii).Dynamics(x_act,u)];
            end
        end

        % Enviromental dynamics
        if isempty(obj.EnviromentModel.TurbulenceModel)
            x_envdot = zeros(0,1);
        else
            x_envdot = obj.EnviromentModel.TurbulenceModel.Dynamics(t,x_env,x_rb);
        end

        % Observer dynamics
        if isempty(obj.Observer)
            x_obsvdot = zeros(0,1);
        else
            x_obsvdot = obj.Observer.Dynamics(t,x_obsv,y);
        end
        
        % Controller dynamics
        if isempty(obj.ControlLaw)
            x_ctrldot = zeros(0,1);
        else
            x_ctrldot = obj.ControlLaw.Dynamics(t,x_ctrl,x_rb,x_obsv,y);
        end
        
        % Assemble state derivative
        xdot = [x_rbdot;x_usdot;x_actdot;x_envdot;x_obsvdot;x_ctrldot];

    end

end % public methods

end