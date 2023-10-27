classdef GNSS_NED < Sensor
    %GNSS_NED Summary of this class goes here
    %

    properties
        % No additional
    end

    methods
        function obj = GNSS_NED(Name)
            %GNSS_NED Construct an instance of this class
            obj.Name = Name;
        end % GNSS_NED

        function y = Model(obj,x,~,~,~)
            %Model Measurement model for a GNSS receiver with measurements
            %in a local NED reference frame. The output, y,
            % is y = [NED;v_NED].

            % Euler angles and rotation matrices
            phi = x(4,1);
            theta = x(5,1);
            psi = x(6,1);
            R1 = [1,0,0;0,cos(phi),-sin(phi);0,sin(phi),cos(phi)];
            R2 = [cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)];
            R3 = [cos(psi),-sin(psi),0;sin(psi),cos(psi),0;0,0,1];
            R_IB = R3*R2*R1;

            % Output
            y = [x(1:3,1);R_IB*x(7:9,1)];

        end % Model

        function [dhdx,dhdu,dhdp] = Jacobian(obj,x,u,Parameters,~)
            %Jacobian Measurement model Jacobians with respect to the state
            % vector, x = [NED; Theta; vb; omega; eta], and the parameter
            % vector, p. Here, eta are any additional states used in the
            % aerodynamic model.

            % Get necessary dimensions.
            nx = size(x,1);
            nu = size(u,1);
            if isobject(Parameters)
                np = size(Parameters.Mean,1);
            else
                np = size(Parameters,1);
            end

            % Euler angles and rotation matrices
            phi = x(4,1);
            theta = x(5,1);
            psi = x(6,1);
            R1 = [1,0,0;0,cos(phi),-sin(phi);0,sin(phi),cos(phi)];
            R2 = [cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)];
            R3 = [cos(psi),-sin(psi),0;sin(psi),cos(psi),0;0,0,1];
            R_IB = R3*R2*R1;
            
            % Compute Jacobians.
            dhdx = zeros(6,nx);
            dhdu = zeros(6,nu);
            dhdp = zeros(6,np);
            dhdx(1:3,1:3) = eye(3);
            dhdx(4:6,7:9) = R_IB;
            dhdx(4,4) = x(8,1)*(sin(phi)*sin(psi)+cos(phi)*cos(psi)*sin(theta))+x(9,1)*(cos(phi)*sin(psi)-cos(psi)*sin(phi)*sin(theta));
            dhdx(4,5) = x(9,1)*cos(phi)*cos(psi)*cos(theta)-x(7,1)*cos(psi)*sin(theta)+x(8,1)*cos(psi)*cos(theta)*sin(phi);
            dhdx(4,6) = x(9,1)*(cos(psi)*sin(phi)-cos(phi)*sin(psi)*sin(theta))-x(8,1)*(cos(phi)*cos(psi)+sin(phi)*sin(psi)*sin(theta))-x(7,1)*cos(theta)*sin(psi);
            dhdx(5,4) = -x(8,1)*(cos(psi)*sin(phi)-cos(phi)*sin(psi)*sin(theta))-x(9,1)*(cos(phi)*cos(psi)+sin(phi)*sin(psi)*sin(theta));
            dhdx(5,5) = x(9,1)*cos(phi)*cos(theta)*sin(psi)-x(7,1)*sin(psi)*sin(theta)+x(8,1)*cos(theta)*sin(phi)*sin(psi);
            dhdx(5,6) = x(9,1)*(sin(phi)*sin(psi)+cos(phi)*cos(psi)*sin(theta))-x(8,1)*(cos(phi)*sin(psi)-cos(psi)*sin(phi)*sin(theta))+x(7,1)*cos(psi)*cos(theta);
            dhdx(6,4) = x(8,1)*cos(phi)*cos(theta)-x(9,1)*cos(theta)*sin(phi);
            dhdx(6,5) = -x(7,1)*cos(theta)-x(9,1)*cos(phi)*sin(theta)-x(8,1)*sin(phi)*sin(theta);

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
                newObj = GNSS_NED(s.Name); 
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