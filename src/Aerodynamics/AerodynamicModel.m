classdef (Abstract) AerodynamicModel < handle & dynamicprops
    %AerodynamicModel Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Name
        Aircraft
        ParameterEstimates
        NumberOfParameters
        NumberOfStates
    end % properties
    
    methods(Abstract)

        % Aerodynamic model structure that returns the bofy frame forces
        % and moments given necessary air-relative aircraft states, any
        % additional states, eta, used in this model, actuator
        % states/inputs delta, aerodynamic parameters, and enviromental
        % constants.
        [Force,Moment] = Model(obj,eta,vb_r,omega_r,delta,Parameters,Constants)

        % Dynamics of additional states, eta.
        detadt = AugmentedDynamics(obj,eta,vb_r,omega_r,delta,Parameters,Constants)

        % Aerodynamic model Jacobians. dfdx is the Jacobian of f with
        % respect to the states, x = [vb;omega;eta;w;omegaw], where eta are
        % any additional states used in this model, w is the wind vector,
        % and omegawis the anguar velocity of the wind as seen by the
        % aircraft. The implementation of this method should use the
        % Aircraft class public method "RigidBodyJacobian". dfdu is the
        % Jacobian of f with respect to the aircraft inputs, u. dfdw is the
        % Jacobian of f with respect to process noise/disturbance, w. dfdp
        % is the Jacobian of f with respect to the aerodynamic parameter
        % vector, p. If any of these are not (yet) required, they may be
        % set as empty arrays.
        [dfdx,dfdu,dfdw,dfdp] = Jacobian(obj,eta,vb_r,omega_r,delta,Parameters,Constants)

        obj = loadobj(s)

    end % Abstract methods

    methods
        
        function obj = AerodynamicModel(Aircraft,Name)
            %AerodynamicModel Construct an instance of this class

            % Required model properties
            obj.Name = Name;
            obj.Aircraft = Aircraft;

            % add to Aircraft
            Aircraft.AddAerodynamicModel(obj,Name);

        end % AerodynamicModel

        function [f,A,D] = Dynamics(obj,t,x,u,w,dervflag,Constants,Parameters)
            %Dynamics Evaluate the dynamic equations of the aircraft for
            % which this aerodynamic model was developed.
            %
            % The states of this system are
            %   x = [NED; Theta; vb; omega; eta]
            % where eta are any additional dynamics used in this model. 
            % The input vector is denoted "u".
            %
            % "Parameters" may be either a ParameterEstimates
            % object or an np x 1 array containing parameter estimates.
            %
            % "Constants" is a struct of which one field must be the
            % gravitational acceleration, g. 

            % If w is empty, set it to zeros.
            nv = 3;
            if isempty(w)
                w = zeros(nv,1);
            end

            % Add in the process noise of wind and angular wind.
            x(7:9,1) = x(7:9,1) + w;
     
            % Compute the modeled forces and moments.
            [F,M] = obj.Model(x,u,Parameters,Constants);
          
            % Evaluate rigid body dynamics.
            dxdt = obj.Aircraft.RigidBodyDynamics(x,F,M,Constants.g);  % !!!! OLD !!!!

            % If there are additional states, their time derivative is
            % computed using the implementation of AugmentedDynamics.
            detadt = obj.AugmentedDynamics(x,u,Parameters,Constants);
            
            % Concatenate state derivatives.
            f = [dxdt; detadt];

            % If the Jacobians are not needed, return empty arrays.
            if dervflag == 0
                A = [];
                D = [];
                return
            end

            % Compute the Jacobians of the aircraft dynamics.
            [A,~,D,~] = obj.Jacobian(x,u,Parameters,Constants);
            
        end % Dynamics

        function [f,A,D] = ParameterDynamics(obj,t,xa,u,vtil,dervflag,Constants)
            %ParameterDynamics Evaluate the dynamic equations of the
            % aircraft augmented with zeros for the dynamics of the
            % aerodynamic parameters for use in a Kalman filter approach to
            % parameter estimation. This function is in the necessary form
            % for use with the Estimation-Tools functions found at
            %       https://github.com/jwhgit01/Estimation-Tools.git
            % These dynamics are in the form
            %
            %                   dxa/dt = f(t,xa,u,vtil)           (1)
            %
            % where vtil(t) is a continuous-time random process.
            %
            % Inputs:
            %
            %   t           The current time in seconds.
            %
            %   xa          The nx x 1 state vector at time t,
            %               xa = [NED; Theta; vb; omega; eta; p]
            %
            %   u           The nu x 1 control vector at time t.
            %
            %   vtil        The nv x 1 process noise vector at time t.
            %
            %   dervflag    A flag that tells whether (dervflag = 1) or not
            %               (dervflag = 0) the partial derivatives df/dx
            %               and df/dvtil must be calculated. If
            %               dervflag = 0, then these outputs will be empty
            %               arrays.
            %
            %   Constants   A structure array of system constants needed to
            %               compute the dynamics. One field must be "g"
            %               containing the gravitational acceleration.
            %  
            % Outputs:
            %
            %   f           The time derivative of x at time t from Eq.(1).
            %
            %   A           The partial derivative of f with respect to x.
            %               This is a Jacobian matrix. It is evaluated and
            %               output only if dervflag = 1.  Otherwise, an
            %               empty array is output.
            %
            %   D           The partial derivative of f with respect to
            %               vtil. This is a Jacobian matrix. It is
            %               evaluated and output only if dervflag = 1.
            %               Otherwise, an empty array is output.
            %

            % Get the necessary dimensions from the Constants input.
            np = Constants.np;
            nxa = size(xa,1);
            nx = nxa-np;

            % Determine the parameter dynamics from the Constants input.
            ParameterModel = Constants.ParameterModel;
            p0 = Constants.p0;

            % Parse the augmented state vector.
            x = xa(1:nx,1);
            p = xa(nx+1:nxa);

            % If vtil is empty, set it to zeros.
            % For aircraft, assume wind is the only process disturbance.
            if strcmp(ParameterModel,'Trivial')
                nv = 3;
            else
                nv = 3 + np;
            end
            if isempty(vtil)
                vtil = zeros(nv,1);
            end

            % Add in the process noise for wind.
            x(7:9,1) = x(7:9,1) + vtil(1:3,1);
            
            % Compute the aircraft dynamics.
            dxdt = obj.Dynamics(x,u,p,Constants);
            
            % Augment the dynamics.
            if strcmp(ParameterModel,'Trivial')
                pdot = zeros(np,1);
            elseif strcmp(ParameterModel,'GaussMarkov')
                pdot = -(1/20)*p + vtil(4:nxa,1);
            elseif strcmp(ParameterModel,'RandomWalk')
                pdot = vtil(4:nxa,1);
            elseif strcmp(ParameterModel,'BiasedGaussMarkov')
                pdot = -(1/20)*(p-p0) + vtil(4:nv,1);
            else
                error('Unknown Parameter Dynamics Model')
            end
            f = [dxdt;pdot];

            % If the Jacobians are not needed, return empty arrays.
            if dervflag == 0
                A = [];
                D = [];
                return
            end

            % Compute the Jacobians of the aircraft dynamics.
            [dfdx,~,dfdw,dfdp] = obj.Jacobian(x,u,p,Constants);

            % Construct the augmented system Jacobians.
            A = [dfdx,dfdp;zeros(np,nxa)];
            if strcmp(ParameterModel,'Trivial')
                D = [dfdw;zeros(np,3)];
            else
                D = [dfdw,zeros(nx,np);zeros(np,3),eye(np)];
            end
            
        end % ParameterDynamics

        function AddParameterEstimates(obj,ParamEstObj)
            %AddParameterEstimates Add parameter estimates to an
            % aerodynamic model.
            obj.ParameterEstimates.(ParamEstObj.Name) = ParamEstObj;
        end % AddParameterEstimates

    end % public methods

end % classdef