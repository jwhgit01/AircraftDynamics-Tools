classdef Magnetometer < Sensor
    %Magnetometer Summary of this class goes here
    %
    % Assumptions:
    %   - The accelerometer is mounted aligned with the aircraft body frame
    %   - The magnetometer bias is constant.

    properties
        Bias
    end

    methods
        function obj = Magnetometer(Name)
            %Accelerometer Construct an instance of this class
            obj.Name = Name;
        end % Accelerometer

        function y = Model(obj,x,~,~,Constants)
            %Model Measurement model for a magnetometer. The output, y,
            % is y = R_BI*BE + bias, where BE is the magnetic field of the
            % Earth at the vehicles current position.
            
            % get variables from Constants struct.
            HomePosition = Constants.HomePosition;
            FlightDate = Constants.FlightDate;

            % Euler angles and rotation matrices
            phi = x(4,1);
            theta = x(5,1);
            psi = x(6,1);
            R1 = [1,0,0;0,cos(phi),-sin(phi);0,sin(phi),cos(phi)];
            R2 = [cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)];
            R3 = [cos(psi),-sin(psi),0;sin(psi),cos(psi),0;0,0,1];

            % rotation matrix from inertial to body frame
            R_BI = (R3*R2*R1).';

            % Compute Earth's magnetic field using the IGRF model
            alt_m = HomePosition(1);
            lat_deg = HomePosition(2);
            lon_deg = HomePosition(3);
            year_dec = decyear(FlightDate);
            BE_nT = igrfmagm(alt_m,lat_deg,lon_deg,year_dec).';
            
            % rotate Earth's magnetic field into the body frame
            y = R_BI*BE_nT;

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

            % get variables from Constants struct.
            HomePosition = Constants.HomePosition;
            FlightDate = Constants.FlightDate;

            % Euler angles
            phi = x(4,1);
            theta = x(5,1);
            psi = x(6,1);

            % Compute Earth's magnetic field using the IGRF model
            alt_m = HomePosition(1);
            lat_deg = HomePosition(2);
            lon_deg = HomePosition(3);
            year_dec = decyear(FlightDate);
            BE_nT = igrfmagm(alt_m,lat_deg,lon_deg,year_dec);

            % Get necessary dimensions.
            nx = size(x,1);
            nu = size(u,1);
            if isobject(Parameters)
                np = size(Parameters.Mean,1);
            else
                np = size(Parameters,1);
            end
            
            % Compute Jacobians.
            dhdx = zeros(3,nx);
            dhdu = zeros(3,nu);
            dhdp = zeros(3,np);
            dhdx(1,5) = -BE_nT(3)*cos(theta)-BE_nT(1)*cos(psi)*sin(theta)-BE_nT(2)*sin(psi)*sin(theta);
            dhdx(1,6) = cos(theta)*(BE_nT(2)*cos(psi)-BE_nT(1)*sin(psi));
            dhdx(2,4) = BE_nT(1)*(sin(phi)*sin(psi)+cos(phi)*cos(psi)*sin(theta))-BE_nT(2)*(cos(psi)*sin(phi)-cos(phi)*sin(psi)*sin(theta))+BE_nT(3)*cos(phi)*cos(theta);
            dhdx(2,5) = sin(phi)*(BE_nT(1)*cos(psi)*cos(theta)-BE_nT(3)*sin(theta)+BE_nT(2)*cos(theta)*sin(psi));
            dhdx(2,6) = -BE_nT(1)*(cos(phi)*cos(psi)+sin(phi)*sin(psi)*sin(theta))-BE_nT(2)*(cos(phi)*sin(psi)-cos(psi)*sin(phi)*sin(theta));
            dhdx(3,4) = BE_nT(1)*(cos(phi)*sin(psi)-cos(psi)*sin(phi)*sin(theta))-BE_nT(2)*(cos(phi)*cos(psi)+sin(phi)*sin(psi)*sin(theta))-BE_nT(3)*cos(theta)*sin(phi);
            dhdx(3,5) = cos(phi)*(BE_nT(1)*cos(psi)*cos(theta)-BE_nT(3)*sin(theta)+BE_nT(2)*cos(theta)*sin(psi));
            dhdx(3,6) = BE_nT(1)*(cos(psi)*sin(phi)-cos(phi)*sin(psi)*sin(theta))+BE_nT(2)*(sin(phi)*sin(psi)+cos(phi)*cos(psi)*sin(theta));

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
                newObj = Magnetometer(s.Name); 
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