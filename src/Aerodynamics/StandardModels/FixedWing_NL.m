classdef FixedWing_NL < AerodynamicModel
    %FixedWing_NL

    properties
        NumProps
        PropEfficientcy
    end

    methods
    
        function obj = FixedWing_NL(Aircraft,Name)
            %FixedWing_NL Construct an instance of this class
            %   Call the AerodynamicModel contstructor.
            obj@AerodynamicModel(Aircraft,Name);

            % Set the number of dynamic states and parameters
            obj.NumberOfStates = 0;
            obj.NumberOfParameters = 36;

            % check to see if the CSV already exists
            AircraftName = obj.Aircraft.Name;
            AircraftPath = obj.Aircraft.WorkspaceDir;
            if 0 ~= exist([AircraftPath '/AircraftData/' AircraftName '_Parameters_FixedWing_NL.csv'],'file')
                warning(['A FixedWing_NL model already exists for ' AircraftName]);
                return
            end
            
            % Copy the template CSV
            temp = which('AircraftDynamcis-Tools/lib/Parameters_FixedWing_NL.csv');
            [libpath,~,~] = fileparts(temp);
            copyfile([libpath filesep 'Parameters_FixedWing_NL.csv'],[AircraftPath '/AircraftData/' AircraftName '_Parameters_FixedWing_NL.csv']);

        end % FixedWing_NL_QS

        function [F,M] = Model(obj,~,vb_r,omega_r,delta,Parameters,Constants)
            %Model

            % Parameters
            b = obj.Aircraft.Span;
            c = obj.Aircraft.Chord;
            S = obj.Aircraft.WingArea;

            % Constants
            rho = Constants.rho;

            % Airspeed and aerodyamic angles
            Vr = sqrt(vb_r(1,1)^2+vb_r(2,1)^2+vb_r(3,1)^2);
            alpha = atan2(vb_r(3,1),vb_r(1,1));
            beta = asin(vb_r(2,1)/Vr);

            % Non-dimensionalized angular rates 
            phat = omega_r(1,1)*b/(2*Vr);
            qhat = omega_r(2,1)*c/(2*Vr);
            rhat = omega_r(3,1)*b/(2*Vr);

            % Actuator states
            da = delta(1,1); % rad
            de = delta(2,1); % rad
            dr = delta(3,1); % rad
            Omega = delta(4,1); % rad/s

            % Parameters
            if isobject(Parameters)
                C = Parameters.Mean;
            else
                C = Parameters;
            end
            %
            % C_x
            %
            C_x_0 = C(1);
            C_x_alpha = C(2);
            C_x_q = C(3);
            C_x_de = C(4);
            C_x_alpha2 = C(5);
            %
            % C_y
            %
            C_y_beta = C(6);
            C_y_p = C(7);
            C_y_r = C(8);
            C_y_da = C(9);
            C_y_dr = C(10);
            C_y_beta3 = C(11);
            %
            % C_z
            %
            C_z_0 = C(12);
            C_z_alpha = C(13);
            C_z_q = C(14);
            C_z_de = C(15);
            C_z_alpha2 = C(16);
            %
            % C_l
            %
            C_l_beta = C(17);
            C_l_p = C(18);
            C_l_r = C(19);
            C_l_da = C(20);
            C_l_dr = C(21);
            C_l_beta3 = C(22);
            %
            % C_m
            %
            C_m_0 = C(23);
            C_m_alpha = C(24);
            C_m_q = C(25);
            C_m_de = C(26);
            C_m_alpha2 = C(27);
            %
            % C_n
            %
            C_n_beta = C(28);
            C_n_p = C(29);
            C_n_r = C(30);
            C_n_da = C(31);
            C_n_dr = C(32);
            C_n_beta3 = C(33);
            %
            % C_T
            %
            C_T_0 = C(34);
            C_T_J = C(35);
            C_T_J2 = C(36);

            % Non-dimensionalized forces and moments
            Cx = C_x_0 + C_x_alpha*alpha + C_x_q*qhat + C_x_de*de ...
                 + C_x_alpha2*alpha^2;
            Cy = C_y_beta*beta + C_y_p*phat + C_y_r*rhat + C_y_da*da ...
                 + C_y_dr*dr + C_y_beta3*beta^3;
            Cz = C_z_0 + C_z_alpha*alpha + C_z_q*qhat + C_z_de*de ...
                 + C_z_alpha2*alpha^2;
            Cl = C_l_beta*beta + C_l_p*phat + C_l_r*rhat + C_l_da*da ...
                 + C_l_dr*dr + C_l_beta3*beta^3;
            Cm = C_m_0 + C_m_alpha*alpha + C_m_q*qhat + C_m_de*de ...
                 + C_m_alpha2*alpha^2;
            Cn = C_n_beta*beta + C_n_p*phat + C_n_r*rhat + C_n_da*da ...
                 + C_n_dr*dr + C_n_beta3*beta^3;

            % Dynamic pressure
            qbar = 0.5*rho*Vr^2;

            % Dimenisonalized aerodynamic forces and moments
            Cxyz = [Cx;Cy;Cz];
            Clmn = [Cl;Cm;Cn];
            Fa = qbar*S*Cxyz;
            Ma = qbar*S*diag([b,c,b])*Clmn;

            % Thrust forces
            D = obj.Aircraft.PropDiameter;
            if Omega < 1
                T = 0;
            else
                J = Vr/(Omega*D);
                CT = C_T_0 + C_T_J*J + C_T_J2*J^2;
                T = obj.PropEfficientcy*obj.NumProps*rho*Omega^2*D^4*CT;
            end

            % Total forces and moments
            F = Fa + [T;0;0];
            M = Ma;
            
        end % Model

        function detadt = AugmentedDynamics(obj,~,~,~,~,~,~)
            detadt = zeros(0,1);
        end

        function [dfdv,dfdomega,dfddelta,dfdp,dfdeta] = Jacobian(obj,~,vb_r,omega_r,delta,Parameters,Constants)
            %Jacobian

            % Get constants and parameters
            b = obj.Aircraft.Span;
            c = obj.Aircraft.Chord;
            S = obj.Aircraft.WingArea;
            D = obj.Aircraft.PropDiameter;
            numprops = obj.NumProps;
            propefficiency = obj.PropEfficientcy;
            rho = Constants.rho;
            if isobject(Parameters)
                C = Parameters.Mean;
            else
                C = Parameters;
            end

            % States and imputs
            u = vb_r(1,1);
            v = vb_r(2,1);
            w = vb_r(3,1);
            p = omega_r(1,1);
            q = omega_r(2,1);
            r = omega_r(3,1);
            da = delta(1,1);
            de = delta(2,1);
            dr = delta(3,1);
            Omega = delta(4,1);

            % The following were computed on 3/22/2023 using symbolic
            % computation in matlab:

            dfdv = [D^4*Omega^2*numprops*propefficiency*rho*((2*C(36)*u)/(D^2*Omega^2)+(C(35)*u)/(D*Omega*(u^2+v^2+w^2)^(1/2)))+S*rho*u*(C(1)+C(4)*de+C(5)*atan2(w,u)^2+C(2)*atan2(w,u)+(C(3)*c*q)/(2*(u^2+v^2+w^2)^(1/2)))-(S*rho*((C(2)*w)/(u^2+w^2)+(2*C(5)*w*atan2(w,u))/(u^2+w^2)+(C(3)*c*q*u)/(2*(u^2+v^2+w^2)^(3/2)))*(u^2+v^2+w^2))/2,C(1)*S*rho*v+C(4)*S*de*rho*v+C(5)*S*rho*v*atan2(w,u)^2+C(2)*S*rho*v*atan2(w,u)+2*C(36)*D^2*numprops*propefficiency*rho*v+(C(3)*S*c*q*rho*v)/(4*(u^2+v^2+w^2)^(1/2))+(C(35)*D^3*Omega*numprops*propefficiency*rho*v)/(u^2+v^2+w^2)^(1/2),(S*rho*((C(2)*u)/(u^2+w^2)+(2*C(5)*u*atan2(w,u))/(u^2+w^2)-(C(3)*c*q*w)/(2*(u^2+v^2+w^2)^(3/2)))*(u^2+v^2+w^2))/2+S*rho*w*(C(1)+C(4)*de+C(5)*atan2(w,u)^2+C(2)*atan2(w,u)+(C(3)*c*q)/(2*(u^2+v^2+w^2)^(1/2)))+D^4*Omega^2*numprops*propefficiency*rho*((2*C(36)*w)/(D^2*Omega^2)+(C(35)*w)/(D*Omega*(u^2+v^2+w^2)^(1/2)));
                C(9)*S*da*rho*u+C(10)*S*dr*rho*u+C(6)*S*rho*u*asin(v/(u^2+v^2+w^2)^(1/2))+C(11)*S*rho*u*asin(v/(u^2+v^2+w^2)^(1/2))^3-(C(6)*S*rho*u*v)/(2*(u^2+w^2)^(1/2))-(3*C(11)*S*rho*u*v*asin(v/(u^2+v^2+w^2)^(1/2))^2)/(2*(u^2+w^2)^(1/2))+(C(7)*S*b*p*rho*u)/(4*(u^2+v^2+w^2)^(1/2))+(C(8)*S*b*r*rho*u)/(4*(u^2+v^2+w^2)^(1/2)),C(9)*S*da*rho*v+C(10)*S*dr*rho*v+C(6)*S*rho*v*asin(v/(u^2+v^2+w^2)^(1/2))+(C(6)*S*rho*u^2)/(2*(u^2+w^2)^(1/2))+(C(6)*S*rho*w^2)/(2*(u^2+w^2)^(1/2))+C(11)*S*rho*v*asin(v/(u^2+v^2+w^2)^(1/2))^3+(3*C(11)*S*rho*u^2*asin(v/(u^2+v^2+w^2)^(1/2))^2)/(2*(u^2+w^2)^(1/2))+(3*C(11)*S*rho*w^2*asin(v/(u^2+v^2+w^2)^(1/2))^2)/(2*(u^2+w^2)^(1/2))+(C(7)*S*b*p*rho*v)/(4*(u^2+v^2+w^2)^(1/2))+(C(8)*S*b*r*rho*v)/(4*(u^2+v^2+w^2)^(1/2)),C(9)*S*da*rho*w+C(10)*S*dr*rho*w+C(6)*S*rho*w*asin(v/(u^2+v^2+w^2)^(1/2))+C(11)*S*rho*w*asin(v/(u^2+v^2+w^2)^(1/2))^3-(C(6)*S*rho*v*w)/(2*(u^2+w^2)^(1/2))-(3*C(11)*S*rho*v*w*asin(v/(u^2+v^2+w^2)^(1/2))^2)/(2*(u^2+w^2)^(1/2))+(C(7)*S*b*p*rho*w)/(4*(u^2+v^2+w^2)^(1/2))+(C(8)*S*b*r*rho*w)/(4*(u^2+v^2+w^2)^(1/2));
                -(S*rho*((C(13)*w)/(u^2+w^2)+(2*C(16)*w*atan2(w,u))/(u^2+w^2)+(C(14)*c*q*u)/(2*(u^2+v^2+w^2)^(3/2)))*(u^2+v^2+w^2))/2+S*rho*u*(C(12)+C(15)*de+C(16)*atan2(w,u)^2+C(13)*atan2(w,u)+(C(14)*c*q)/(2*(u^2+v^2+w^2)^(1/2))),C(12)*S*rho*v+C(15)*S*de*rho*v+C(16)*S*rho*v*atan2(w,u)^2+C(13)*S*rho*v*atan2(w,u)+(C(14)*S*c*q*rho*v)/(4*(u^2+v^2+w^2)^(1/2)),(S*rho*((C(13)*u)/(u^2+w^2)+(2*C(16)*u*atan2(w,u))/(u^2+w^2)-(C(14)*c*q*w)/(2*(u^2+v^2+w^2)^(3/2)))*(u^2+v^2+w^2))/2+S*rho*w*(C(12)+C(15)*de+C(16)*atan2(w,u)^2+C(13)*atan2(w,u)+(C(14)*c*q)/(2*(u^2+v^2+w^2)^(1/2)));
                C(22)*S*b*rho*u*asin(v/(u^2+v^2+w^2)^(1/2))^3+C(20)*S*b*da*rho*u+C(21)*S*b*dr*rho*u+C(17)*S*b*rho*u*asin(v/(u^2+v^2+w^2)^(1/2))-(C(17)*S*b*rho*u*v)/(2*(u^2+w^2)^(1/2))+(C(18)*S*b^2*p*rho*u)/(4*(u^2+v^2+w^2)^(1/2))+(C(19)*S*b^2*r*rho*u)/(4*(u^2+v^2+w^2)^(1/2))-(3*C(22)*S*b*rho*u*v*asin(v/(u^2+v^2+w^2)^(1/2))^2)/(2*(u^2+w^2)^(1/2)),(S*b*rho*(2*C(17)*u^2+2*C(17)*w^2+6*C(22)*u^2*asin(v/(u^2+v^2+w^2)^(1/2))^2+6*C(22)*w^2*asin(v/(u^2+v^2+w^2)^(1/2))^2-(C(18)*b*p*v*(u^2+w^2)^(1/2))/(u^2+v^2+w^2)^(1/2)-(C(19)*b*r*v*(u^2+w^2)^(1/2))/(u^2+v^2+w^2)^(1/2)))/(4*(u^2+w^2)^(1/2))+S*b*rho*v*(C(20)*da+C(21)*dr+C(17)*asin(v/(u^2+v^2+w^2)^(1/2))+C(22)*asin(v/(u^2+v^2+w^2)^(1/2))^3+(C(18)*b*p)/(2*(u^2+v^2+w^2)^(1/2))+(C(19)*b*r)/(2*(u^2+v^2+w^2)^(1/2))),C(22)*S*b*rho*w*asin(v/(u^2+v^2+w^2)^(1/2))^3+C(20)*S*b*da*rho*w+C(21)*S*b*dr*rho*w+C(17)*S*b*rho*w*asin(v/(u^2+v^2+w^2)^(1/2))-(C(17)*S*b*rho*v*w)/(2*(u^2+w^2)^(1/2))+(C(18)*S*b^2*p*rho*w)/(4*(u^2+v^2+w^2)^(1/2))+(C(19)*S*b^2*r*rho*w)/(4*(u^2+v^2+w^2)^(1/2))-(3*C(22)*S*b*rho*v*w*asin(v/(u^2+v^2+w^2)^(1/2))^2)/(2*(u^2+w^2)^(1/2));
                -(S*c*rho*((C(24)*w)/(u^2+w^2)+(2*C(27)*w*atan2(w,u))/(u^2+w^2)+(C(25)*c*q*u)/(2*(u^2+v^2+w^2)^(3/2)))*(u^2+v^2+w^2))/2+S*c*rho*u*(C(23)+C(26)*de+C(27)*atan2(w,u)^2+C(24)*atan2(w,u)+(C(25)*c*q)/(2*(u^2+v^2+w^2)^(1/2))),C(23)*S*c*rho*v+C(24)*S*c*rho*v*atan2(w,u)+C(26)*S*c*de*rho*v+C(27)*S*c*rho*v*atan2(w,u)^2+(C(25)*S*c^2*q*rho*v)/(4*(u^2+v^2+w^2)^(1/2)),(S*c*rho*((C(24)*u)/(u^2+w^2)+(2*C(27)*u*atan2(w,u))/(u^2+w^2)-(C(25)*c*q*w)/(2*(u^2+v^2+w^2)^(3/2)))*(u^2+v^2+w^2))/2+S*c*rho*w*(C(23)+C(26)*de+C(27)*atan2(w,u)^2+C(24)*atan2(w,u)+(C(25)*c*q)/(2*(u^2+v^2+w^2)^(1/2)));
                C(33)*S*b*rho*u*asin(v/(u^2+v^2+w^2)^(1/2))^3+C(31)*S*b*da*rho*u+C(32)*S*b*dr*rho*u+C(28)*S*b*rho*u*asin(v/(u^2+v^2+w^2)^(1/2))-(C(28)*S*b*rho*u*v)/(2*(u^2+w^2)^(1/2))+(C(29)*S*b^2*p*rho*u)/(4*(u^2+v^2+w^2)^(1/2))+(C(30)*S*b^2*r*rho*u)/(4*(u^2+v^2+w^2)^(1/2))-(3*C(33)*S*b*rho*u*v*asin(v/(u^2+v^2+w^2)^(1/2))^2)/(2*(u^2+w^2)^(1/2)),(S*b*rho*(2*C(28)*u^2+2*C(28)*w^2+6*C(33)*u^2*asin(v/(u^2+v^2+w^2)^(1/2))^2+6*C(33)*w^2*asin(v/(u^2+v^2+w^2)^(1/2))^2-(C(29)*b*p*v*(u^2+w^2)^(1/2))/(u^2+v^2+w^2)^(1/2)-(C(30)*b*r*v*(u^2+w^2)^(1/2))/(u^2+v^2+w^2)^(1/2)))/(4*(u^2+w^2)^(1/2))+S*b*rho*v*(C(31)*da+C(32)*dr+C(28)*asin(v/(u^2+v^2+w^2)^(1/2))+C(33)*asin(v/(u^2+v^2+w^2)^(1/2))^3+(C(29)*b*p)/(2*(u^2+v^2+w^2)^(1/2))+(C(30)*b*r)/(2*(u^2+v^2+w^2)^(1/2))),C(33)*S*b*rho*w*asin(v/(u^2+v^2+w^2)^(1/2))^3+C(31)*S*b*da*rho*w+C(32)*S*b*dr*rho*w+C(28)*S*b*rho*w*asin(v/(u^2+v^2+w^2)^(1/2))-(C(28)*S*b*rho*v*w)/(2*(u^2+w^2)^(1/2))+(C(29)*S*b^2*p*rho*w)/(4*(u^2+v^2+w^2)^(1/2))+(C(30)*S*b^2*r*rho*w)/(4*(u^2+v^2+w^2)^(1/2))-(3*C(33)*S*b*rho*v*w*asin(v/(u^2+v^2+w^2)^(1/2))^2)/(2*(u^2+w^2)^(1/2))];
                
                
            dfdomega = [0,(C(3)*S*c*rho*(u^2+v^2+w^2)^(1/2))/4,0;
                (C(7)*S*b*rho*(u^2+v^2+w^2)^(1/2))/4,0,(C(8)*S*b*rho*(u^2+v^2+w^2)^(1/2))/4;
                0,(C(14)*S*c*rho*(u^2+v^2+w^2)^(1/2))/4,0;
                (C(18)*S*b^2*rho*(u^2+v^2+w^2)^(1/2))/4,0,(C(19)*S*b^2*rho*(u^2+v^2+w^2)^(1/2))/4;
                0,(C(25)*S*c^2*rho*(u^2+v^2+w^2)^(1/2))/4,0;
                (C(29)*S*b^2*rho*(u^2+v^2+w^2)^(1/2))/4,0,(C(30)*S*b^2*rho*(u^2+v^2+w^2)^(1/2))/4];
                
                
            dfddelta = [0,(C(4)*S*rho*(u^2+v^2+w^2))/2,0,D^3*numprops*propefficiency*rho*(C(35)*(u^2+v^2+w^2)^(1/2)+2*C(34)*D*Omega);
                (C(9)*S*rho*(u^2+v^2+w^2))/2,0,(C(10)*S*rho*(u^2+v^2+w^2))/2,0;
                0,(C(15)*S*rho*(u^2+v^2+w^2))/2,0,0;
                (C(20)*S*b*rho*(u^2+v^2+w^2))/2,0,(C(21)*S*b*rho*(u^2+v^2+w^2))/2,0;
                0,(C(26)*S*c*rho*(u^2+v^2+w^2))/2,0,0;
                (C(31)*S*b*rho*(u^2+v^2+w^2))/2,0,(C(32)*S*b*rho*(u^2+v^2+w^2))/2,0];
                
                
            dfdp = [(S*rho*(u^2+v^2+w^2))/2,(S*rho*atan2(w,u)*(u^2+v^2+w^2))/2,(S*c*q*rho*(u^2+v^2+w^2)^(1/2))/4,(S*de*rho*(u^2+v^2+w^2))/2,(S*rho*atan2(w,u)^2*(u^2+v^2+w^2))/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,D^4*Omega^2*numprops*propefficiency*rho,D^3*Omega*numprops*propefficiency*rho*(u^2+v^2+w^2)^(1/2),D^2*numprops*propefficiency*rho*(u^2+v^2+w^2);
                0,0,0,0,0,(S*rho*asin(v/(u^2+v^2+w^2)^(1/2))*(u^2+v^2+w^2))/2,(S*b*p*rho*(u^2+v^2+w^2)^(1/2))/4,(S*b*r*rho*(u^2+v^2+w^2)^(1/2))/4,(S*da*rho*(u^2+v^2+w^2))/2,(S*dr*rho*(u^2+v^2+w^2))/2,(S*rho*asin(v/(u^2+v^2+w^2)^(1/2))^3*(u^2+v^2+w^2))/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                0,0,0,0,0,0,0,0,0,0,0,(S*rho*(u^2+v^2+w^2))/2,(S*rho*atan2(w,u)*(u^2+v^2+w^2))/2,(S*c*q*rho*(u^2+v^2+w^2)^(1/2))/4,(S*de*rho*(u^2+v^2+w^2))/2,(S*rho*atan2(w,u)^2*(u^2+v^2+w^2))/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(S*b*rho*asin(v/(u^2+v^2+w^2)^(1/2))*(u^2+v^2+w^2))/2,(S*b^2*p*rho*(u^2+v^2+w^2)^(1/2))/4,(S*b^2*r*rho*(u^2+v^2+w^2)^(1/2))/4,(S*b*da*rho*(u^2+v^2+w^2))/2,(S*b*dr*rho*(u^2+v^2+w^2))/2,(S*b*rho*asin(v/(u^2+v^2+w^2)^(1/2))^3*(u^2+v^2+w^2))/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(S*c*rho*(u^2+v^2+w^2))/2,(S*c*rho*atan2(w,u)*(u^2+v^2+w^2))/2,(S*c^2*q*rho*(u^2+v^2+w^2)^(1/2))/4,(S*c*de*rho*(u^2+v^2+w^2))/2,(S*c*rho*atan2(w,u)^2*(u^2+v^2+w^2))/2,0,0,0,0,0,0,0,0,0;
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(S*b*rho*asin(v/(u^2+v^2+w^2)^(1/2))*(u^2+v^2+w^2))/2,(S*b^2*p*rho*(u^2+v^2+w^2)^(1/2))/4,(S*b^2*r*rho*(u^2+v^2+w^2)^(1/2))/4,(S*b*da*rho*(u^2+v^2+w^2))/2,(S*b*dr*rho*(u^2+v^2+w^2))/2,(S*b*rho*asin(v/(u^2+v^2+w^2)^(1/2))^3*(u^2+v^2+w^2))/2,0,0,0];

            dfdeta = zeros(6,0);
        
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
                newObj = FixedWing_NL(s.Aircraft,s.Name); 
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