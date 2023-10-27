classdef Accelerometer < Sensor
    %Accelerometer Summary of this class goes here
    %
    % Assumptions:
    %   - The accelerometer is mounted at the center of gravity aligned
    %     with the aircraft body frame.
    %   - The accelerometer bias is constant.

    properties
        AerodynamicModel
        Bias
    end

    methods
        function obj = Accelerometer(Name,AerodynamicModel)
            %Accelerometer Construct an instance of this class
            obj.Name = Name;
            obj.AerodynamicModel = AerodynamicModel;
        end % Accelerometer

        function y = Model(obj,x,u,Parameters,Constants)
            %Model Measurement model for an accelerometer. The output, y,
            % is y = [ax; ay; az] + bias. Note this model assumes the
            % accelerometer frame is oriented with the body frame of the
            % aircraft.
            
            % get variables from Constants struct.
            g = Constants.g;

            % aircraft inertial properties
            m = obj.AerodynamicModel.Aircraft.Mass;

            % Euler angles and rotation matrices
            phi = x(4,1);
            theta = x(5,1);
            psi = x(6,1);
            R1 = [1,0,0;0,cos(phi),-sin(phi);0,sin(phi),cos(phi)];
            R2 = [cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)];
            R3 = [cos(psi),-sin(psi),0;sin(psi),cos(psi),0;0,0,1];

            % rotation matrix from inertial to body frame
            R_BI = (R3*R2*R1).';
            
            % Modeled aerodynamic forces
            [F,~] = obj.AerodynamicModel.Model(x,u,Parameters,Constants);

            % total body acceleration
            e3 = [0;0;1];
            y = g*R_BI*e3 + F/m;

            % If there is bias, add it.
            if ~isempty(obj.Bias)
                y = y + obj.Bias;
            end

        end % Model

        function [dhdx,dhdu,dhdp] = Jacobian(obj,x,u,Parameters,Constants)
            %Jacobian Measurement model Jacobians with respect to the state
            % vector, x = [NED; Theta; vb; omega; eta], and the parameter
            % vector, p. Here, eta are any additional states used in the
            % aerodynamic model.
            
            % Compute Jacobians of aircraft dynamics
            [dfdx,dfdu,~,dfdp] = obj.AerodynamicModel.Jacobian(x,u,Parameters,Constants);

            % Get necessary dimensions
            nx = size(dfdx,2);

            % Compute Jacobian of cross product term in rigid body dynamics
            vb = x(6:9,1);
            omega = x(10:12,1);
            hat = @(a) [0,-a(3),a(2);a(3),0,-a(1);-a(2),a(1),0];
            dcpdx = [zeros(3,6), -hat(omega), hat(vb), zeros(3,nx-12)];

            % Subtract out cross product term from aircraft dynamics Jacobian.
            dhdx = dfdx(7:9,:) - dcpdx;

            % Assemble other Jacobians
            dhdu = dfdu(7:9,:);
            dhdp = dfdp(7:9,:);

        end % Jacobian
      
    end % public methods

    methods(Static)
        
        function obj = loadobj(s)
            %loadobj
            %
            % If any of the saved properties no longer are defined in the
            % class definition, create a dynamic property and give a
            % warning.
            if isstruct(s)
                newObj = Accelerometer(s.Name,s.AerodynamicModel); 
                f = fieldnames(s);
                for ii = 1:length(f)
                    if isprop(newObj,f{ii})
                        newObj.(f{ii}) = s.(f{ii});
                    else
                        addprop(newObj,f{ii});
                        newObj.(f{ii}) = s.(f{ii});
                        mc = metaclass(newObj);
                        warning(['Property ''' f{ii} ''' loaded into object'...
                            'as dynamic property.\nIt should be defined in '...
                            mc.Name '.m or its superclasses\nin order to '...
                            'maintain compatibility.'],[])
                    end
                end
                obj = newObj;
            else
                obj = s;
            end
        end

    end % static methods

end % classdef