classdef FixedWing_NL_US < AerodynamicModel
    %FixedWing_NL_US

    properties
        UnsteadyTransferFunctions
        UnsteadyStateSpace
        NumProps
        PropEfficientcy
    end

    methods
    
        function obj = FixedWing_NL_US(Aircraft,Name)
            %FixedWing_NL_QS Construct an instance of this class
            %   Call the AerodynamicModel contstructor.
            obj@AerodynamicModel(Aircraft,Name);

            % Initialize UnsteadyTransferFunctions struct
            obj.UnsteadyTransferFunctions.X = [];
            obj.UnsteadyTransferFunctions.Y = [];
            obj.UnsteadyTransferFunctions.Z = [];
            obj.UnsteadyTransferFunctions.L = [];
            obj.UnsteadyTransferFunctions.M = [];
            obj.UnsteadyTransferFunctions.N = [];

            % check to see if the CSV already exists
            AircraftName = obj.Aircraft.Name;
            AircraftPath = obj.Aircraft.WorkspaceDir;
            if 0 ~= exist([AircraftPath '/AircraftData/' AircraftName '_Parameters_FixedWing_NL_US.csv'],'file')
                warning(['A FixedWing_NL_US model already exists for ' AircraftName]);
                return
            end
            
            % Copy the template CSV
            temp = which('AircraftDynamcis-Tools/lib/Parameters_FixedWing_NL_US.csv');
            [libpath,~,~] = fileparts(temp);
            copyfile([libpath filesep 'Parameters_FixedWing_NL_QS.csv'],[AircraftPath '/AircraftData/' AircraftName '_Parameters_FixedWing_NL_US.csv']);

        end % FixedWing_NL_QS

        function [F,M] = Model(obj,eta,vb_r,omega_r,delta,Parameters,Constants)
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

            % Propeller advance ratio
            D = obj.Aircraft.PropDiameter;
            J = Vr/(Omega*D);

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
            C_x_alphadot = C(6);
            %
            % C_y
            %
            C_y_beta = C(7);
            C_y_p = C(8);
            C_y_r = C(9);
            C_y_da = C(10);
            C_y_dr = C(11);
            C_y_beta3 = C(12);
            C_y_betadot = C(13);
            %
            % C_z
            %
            C_z_0 = C(14);
            C_z_alpha = C(15);
            C_z_q = C(16);
            C_z_de = C(17);
            C_z_alpha2 = C(18);
            C_z_alphadot = C(19);
            %
            % C_l
            %
            C_l_beta = C(20);
            C_l_p = C(21);
            C_l_r = C(22);
            C_l_da = C(23);
            C_l_dr = C(24);
            C_l_beta3 = C(25);
            C_l_betadot = C(26);
            %
            % C_m
            %
            C_m_0 = C(27);
            C_m_alpha = C(28);
            C_m_q = C(29);
            C_m_de = C(30);
            C_m_alpha2 = C(31);
            C_m_alphadot = C(32);
            %
            % C_n
            %
            C_n_beta = C(33);
            C_n_p = C(34);
            C_n_r = C(35);
            C_n_da = C(36);
            C_n_dr = C(37);
            C_n_beta3 = C(38);
            C_n_betadot = C(39);
            %
            % C_T
            %
            C_T_0 = C(40);
            C_T_J = C(41);
            C_T_J2 = C(42);

            % Unsteady aerodynamics LTI system
            sys = obj.UnsteadyStateSpace;
            u = C_x_alpha

            % Non-dimensionalized forces and moments
            Cx = C_x_0 + C_x_alpha*alpha + C_x_q*qhat + C_x_de*de ...
                 + C_x_alpha2*alpha^2 + C_x_alphadot*alphadot;
            Cy = C_y_beta*beta + C_y_p*phat + C_y_r*rhat + C_y_da*da ...
                 + C_y_dr*dr + C_y_beta3*beta^3 + C_y_betadot*betadot;
            Cz = C_z_0 + C_z_alpha*alpha + C_z_q*qhat + C_z_de*de ...
                 + C_z_alpha2*alpha^2 + C_z_alphadot*alphadot;
            Cl = C_l_beta*beta + C_l_p*phat + C_l_r*rhat + C_l_da*da ...
                 + C_l_dr*dr + C_l_beta3*beta^3 + C_l_betadot*betadot;
            Cm = C_m_0 + C_m_alpha*alpha + C_m_q*qhat + C_m_de*de ...
                 + C_m_alpha2*alpha^2 + C_m_alphadot*alphadot;
            Cn = C_n_beta*beta + C_n_p*phat + C_n_r*rhat + C_n_da*da ...
                 + C_n_dr*dr + C_n_beta3*beta^3 + C_n_betadot*betadot;

            % Dynamic pressure
            qbar = 0.5*rho*Vr^2;

            % Dimenisonalized aerodynamic forces and moments
            Cxyz = [Cx;Cy;Cz];
            Clmn = [Cl;Cm;Cn];
            Fa = qbar*S*Cxyz;
            Ma = qbar*S*diag([b,c,b])*Clmn;

            % Thrust forces
            CT = C_T_0 + C_T_J*J + C_T_J2*J^2;
            T = obj.PropEfficientcy*obj.NumProps*rho*Omega^2*D^4*CT;

            % Total forces and moments
            F = Fa + [T;0;0];
            M = Ma;
            
        end % Model

        function detadt = AugmentedDynamics(obj,eta,vb_r,~,~,~,~)
            tau = obj.DerivativeTimeConstant;
            Vr = norm(vb_r);
            alpha = atan2(vb_r(3,1),vb_r(1,1));
            beta = asin(vb_r(2,1)/Vr);
            alphabeta = [alpha;beta];
            detadt = -(1/tau)*eye(2)*eta + alphabeta;
        end

        function [dfdv,dfdomega,dfddelta,dfdp] = Jacobian(obj,eta,vb_r,omega_r,delta,Parameters,Constants)
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

            % Airspeed and aerodyamic angles
            Vr = sqrt(vb_r(1,1)^2+vb_r(2,1)^2+vb_r(3,1)^2);
            alpha = atan2(vb_r(3,1),vb_r(1,1));
            beta = asin(vb_r(2,1)/Vr);

            % Quasi-steady aerodyamic angles
            tau = obj.DerivativeTimeConstant;
            alphabeta = [alpha;beta];
            alphabetadot = -1/(tau^2)*eye(2)*eta + 1/tau*eye(2)*alphabeta;
            alphadot = alphabetadot(1);
            betadot = alphabetadot(2);

            % The following were computed on 3/15/2023 using symbolic
            % computation in matlab:

            dfdv = [...
                D^4*Omega^2*numprops*propefficiency*rho*((2*C(42)*u)/(D^2*Omega^2)+(C(41)*u)/(D*Omega*(u^2+v^2+w^2)^(1/2)))+S*rho*u*(C(1)+C(6)*alphadot+C(4)*de+C(5)*atan2(w,u)^2+C(2)*atan2(w,u)+(C(3)*c*q)/(2*(u^2+v^2+w^2)^(1/2)))-(S*rho*((C(2)*w)/(u^2+w^2)+(2*C(5)*w*atan2(w,u))/(u^2+w^2)+(C(3)*c*q*u)/(2*(u^2+v^2+w^2)^(3/2)))*(u^2+v^2+w^2))/2, ...
	                C(1)*S*rho*v+C(6)*S*alphadot*rho*v+C(4)*S*de*rho*v+C(5)*S*rho*v*atan2(w,u)^2+C(2)*S*rho*v*atan2(w,u)+2*C(42)*D^2*numprops*propefficiency*rho*v+(C(3)*S*c*q*rho*v)/(4*(u^2+v^2+w^2)^(1/2))+(C(41)*D^3*Omega*numprops*propefficiency*rho*v)/(u^2+v^2+w^2)^(1/2), ...
	                (S*rho*((C(2)*u)/(u^2+w^2)+(2*C(5)*u*atan2(w,u))/(u^2+w^2)-(C(3)*c*q*w)/(2*(u^2+v^2+w^2)^(3/2)))*(u^2+v^2+w^2))/2+S*rho*w*(C(1)+C(6)*alphadot+C(4)*de+C(5)*atan2(w,u)^2+C(2)*atan2(w,u)+(C(3)*c*q)/(2*(u^2+v^2+w^2)^(1/2)))+D^4*Omega^2*numprops*propefficiency*rho*((2*C(42)*w)/(D^2*Omega^2)+(C(41)*w)/(D*Omega*(u^2+v^2+w^2)^(1/2)));
                C(13)*S*betadot*rho*u+C(10)*S*da*rho*u+C(11)*S*dr*rho*u+C(7)*S*rho*u*asin(v/(u^2+v^2+w^2)^(1/2))+C(12)*S*rho*u*asin(v/(u^2+v^2+w^2)^(1/2))^3-(C(7)*S*rho*u*v)/(2*(u^2+w^2)^(1/2))-(3*C(12)*S*rho*u*v*asin(v/(u^2+v^2+w^2)^(1/2))^2)/(2*(u^2+w^2)^(1/2))+(C(8)*S*b*p*rho*u)/(4*(u^2+v^2+w^2)^(1/2))+(C(9)*S*b*r*rho*u)/(4*(u^2+v^2+w^2)^(1/2)), ...
	                S*rho*v*(C(13)*betadot+C(10)*da+C(11)*dr+C(7)*asin(v/(u^2+v^2+w^2)^(1/2))+C(12)*asin(v/(u^2+v^2+w^2)^(1/2))^3+(C(8)*b*p)/(2*(u^2+v^2+w^2)^(1/2))+(C(9)*b*r)/(2*(u^2+v^2+w^2)^(1/2)))+(S*rho*(2*C(7)*u^2+2*C(7)*w^2+6*C(12)*u^2*asin(v/(u^2+v^2+w^2)^(1/2))^2+6*C(12)*w^2*asin(v/(u^2+v^2+w^2)^(1/2))^2-(C(8)*b*p*v*(u^2+w^2)^(1/2))/(u^2+v^2+w^2)^(1/2)-(C(9)*b*r*v*(u^2+w^2)^(1/2))/(u^2+v^2+w^2)^(1/2)))/(4*(u^2+w^2)^(1/2)), ...
	                C(13)*S*betadot*rho*w+C(10)*S*da*rho*w+C(11)*S*dr*rho*w+C(7)*S*rho*w*asin(v/(u^2+v^2+w^2)^(1/2))+C(12)*S*rho*w*asin(v/(u^2+v^2+w^2)^(1/2))^3-(C(7)*S*rho*v*w)/(2*(u^2+w^2)^(1/2))-(3*C(12)*S*rho*v*w*asin(v/(u^2+v^2+w^2)^(1/2))^2)/(2*(u^2+w^2)^(1/2))+(C(8)*S*b*p*rho*w)/(4*(u^2+v^2+w^2)^(1/2))+(C(9)*S*b*r*rho*w)/(4*(u^2+v^2+w^2)^(1/2));
                -(S*rho*((C(15)*w)/(u^2+w^2)+(2*C(18)*w*atan2(w,u))/(u^2+w^2)+(C(16)*c*q*u)/(2*(u^2+v^2+w^2)^(3/2)))*(u^2+v^2+w^2))/2+S*rho*u*(C(14)+C(19)*alphadot+C(17)*de+C(18)*atan2(w,u)^2+C(15)*atan2(w,u)+(C(16)*c*q)/(2*(u^2+v^2+w^2)^(1/2))), ...
	                C(14)*S*rho*v+C(19)*S*alphadot*rho*v+C(17)*S*de*rho*v+C(18)*S*rho*v*atan2(w,u)^2+C(15)*S*rho*v*atan2(w,u)+(C(16)*S*c*q*rho*v)/(4*(u^2+v^2+w^2)^(1/2)), ...
	                (S*rho*((C(15)*u)/(u^2+w^2)+(2*C(18)*u*atan2(w,u))/(u^2+w^2)-(C(16)*c*q*w)/(2*(u^2+v^2+w^2)^(3/2)))*(u^2+v^2+w^2))/2+S*rho*w*(C(14)+C(19)*alphadot+C(17)*de+C(18)*atan2(w,u)^2+C(15)*atan2(w,u)+(C(16)*c*q)/(2*(u^2+v^2+w^2)^(1/2)));
                C(25)*S*b*rho*u*asin(v/(u^2+v^2+w^2)^(1/2))^3+C(26)*S*b*betadot*rho*u+C(23)*S*b*da*rho*u+C(24)*S*b*dr*rho*u+C(20)*S*b*rho*u*asin(v/(u^2+v^2+w^2)^(1/2))-(C(20)*S*b*rho*u*v)/(2*(u^2+w^2)^(1/2))+(C(21)*S*b^2*p*rho*u)/(4*(u^2+v^2+w^2)^(1/2))+(C(22)*S*b^2*r*rho*u)/(4*(u^2+v^2+w^2)^(1/2))-(3*C(25)*S*b*rho*u*v*asin(v/(u^2+v^2+w^2)^(1/2))^2)/(2*(u^2+w^2)^(1/2)), ...
	                (S*b*rho*(2*C(20)*u^2+2*C(20)*w^2+6*C(25)*u^2*asin(v/(u^2+v^2+w^2)^(1/2))^2+6*C(25)*w^2*asin(v/(u^2+v^2+w^2)^(1/2))^2-(C(21)*b*p*v*(u^2+w^2)^(1/2))/(u^2+v^2+w^2)^(1/2)-(C(22)*b*r*v*(u^2+w^2)^(1/2))/(u^2+v^2+w^2)^(1/2)))/(4*(u^2+w^2)^(1/2))+S*b*rho*v*(C(26)*betadot+C(23)*da+C(24)*dr+C(20)*asin(v/(u^2+v^2+w^2)^(1/2))+C(25)*asin(v/(u^2+v^2+w^2)^(1/2))^3+(C(21)*b*p)/(2*(u^2+v^2+w^2)^(1/2))+(C(22)*b*r)/(2*(u^2+v^2+w^2)^(1/2))), ...
	                C(25)*S*b*rho*w*asin(v/(u^2+v^2+w^2)^(1/2))^3+C(26)*S*b*betadot*rho*w+C(23)*S*b*da*rho*w+C(24)*S*b*dr*rho*w+C(20)*S*b*rho*w*asin(v/(u^2+v^2+w^2)^(1/2))-(C(20)*S*b*rho*v*w)/(2*(u^2+w^2)^(1/2))+(C(21)*S*b^2*p*rho*w)/(4*(u^2+v^2+w^2)^(1/2))+(C(22)*S*b^2*r*rho*w)/(4*(u^2+v^2+w^2)^(1/2))-(3*C(25)*S*b*rho*v*w*asin(v/(u^2+v^2+w^2)^(1/2))^2)/(2*(u^2+w^2)^(1/2));
                -(S*c*rho*((C(28)*w)/(u^2+w^2)+(2*C(31)*w*atan2(w,u))/(u^2+w^2)+(C(29)*c*q*u)/(2*(u^2+v^2+w^2)^(3/2)))*(u^2+v^2+w^2))/2+S*c*rho*u*(C(27)+C(32)*alphadot+C(30)*de+C(31)*atan2(w,u)^2+C(28)*atan2(w,u)+(C(29)*c*q)/(2*(u^2+v^2+w^2)^(1/2))), ...
	                S*c*rho*v*(C(27)+C(32)*alphadot+C(30)*de+C(31)*atan2(w,u)^2+C(28)*atan2(w,u)+(C(29)*c*q)/(2*(u^2+v^2+w^2)^(1/2)))-(C(29)*S*c^2*q*rho*v)/(4*(u^2+v^2+w^2)^(1/2)), ...
	                (S*c*rho*((C(28)*u)/(u^2+w^2)+(2*C(31)*u*atan2(w,u))/(u^2+w^2)-(C(29)*c*q*w)/(2*(u^2+v^2+w^2)^(3/2)))*(u^2+v^2+w^2))/2+S*c*rho*w*(C(27)+C(32)*alphadot+C(30)*de+C(31)*atan2(w,u)^2+C(28)*atan2(w,u)+(C(29)*c*q)/(2*(u^2+v^2+w^2)^(1/2)));
                C(38)*S*b*rho*u*asin(v/(u^2+v^2+w^2)^(1/2))^3+C(39)*S*b*betadot*rho*u+C(36)*S*b*da*rho*u+C(37)*S*b*dr*rho*u+C(33)*S*b*rho*u*asin(v/(u^2+v^2+w^2)^(1/2))-(C(33)*S*b*rho*u*v)/(2*(u^2+w^2)^(1/2))+(C(34)*S*b^2*p*rho*u)/(4*(u^2+v^2+w^2)^(1/2))+(C(35)*S*b^2*r*rho*u)/(4*(u^2+v^2+w^2)^(1/2))-(3*C(38)*S*b*rho*u*v*asin(v/(u^2+v^2+w^2)^(1/2))^2)/(2*(u^2+w^2)^(1/2)), ...
	                (S*b*rho*(2*C(33)*u^2+2*C(33)*w^2+6*C(38)*u^2*asin(v/(u^2+v^2+w^2)^(1/2))^2+6*C(38)*w^2*asin(v/(u^2+v^2+w^2)^(1/2))^2-(C(34)*b*p*v*(u^2+w^2)^(1/2))/(u^2+v^2+w^2)^(1/2)-(C(35)*b*r*v*(u^2+w^2)^(1/2))/(u^2+v^2+w^2)^(1/2)))/(4*(u^2+w^2)^(1/2))+S*b*rho*v*(C(39)*betadot+C(36)*da+C(37)*dr+C(33)*asin(v/(u^2+v^2+w^2)^(1/2))+C(38)*asin(v/(u^2+v^2+w^2)^(1/2))^3+(C(34)*b*p)/(2*(u^2+v^2+w^2)^(1/2))+(C(35)*b*r)/(2*(u^2+v^2+w^2)^(1/2))), ...
	                C(38)*S*b*rho*w*asin(v/(u^2+v^2+w^2)^(1/2))^3+C(39)*S*b*betadot*rho*w+C(36)*S*b*da*rho*w+C(37)*S*b*dr*rho*w+C(33)*S*b*rho*w*asin(v/(u^2+v^2+w^2)^(1/2))-(C(33)*S*b*rho*v*w)/(2*(u^2+w^2)^(1/2))+(C(34)*S*b^2*p*rho*w)/(4*(u^2+v^2+w^2)^(1/2))+(C(35)*S*b^2*r*rho*w)/(4*(u^2+v^2+w^2)^(1/2))-(3*C(38)*S*b*rho*v*w*asin(v/(u^2+v^2+w^2)^(1/2))^2)/(2*(u^2+w^2)^(1/2))];

            dfdomega = [0,(C(3)*S*c*rho*(u^2+v^2+w^2)^(1/2))/4,0;
                (C(8)*S*b*rho*(u^2+v^2+w^2)^(1/2))/4,0,(C(9)*S*b*rho*(u^2+v^2+w^2)^(1/2))/4;
                0,(C(16)*S*c*rho*(u^2+v^2+w^2)^(1/2))/4,0;
                (C(21)*S*b^2*rho*(u^2+v^2+w^2)^(1/2))/4,0,(C(22)*S*b^2*rho*(u^2+v^2+w^2)^(1/2))/4;
                0,(C(29)*S*c^2*rho*(u^2+v^2+w^2)^(1/2))/4,0;
                (C(34)*S*b^2*rho*(u^2+v^2+w^2)^(1/2))/4,0,(C(35)*S*b^2*rho*(u^2+v^2+w^2)^(1/2))/4];
            
            dfddelta = [0,(C(4)*S*rho*(u^2+v^2+w^2))/2,0,D^3*numprops*propefficiency*rho*(C(41)*(u^2+v^2+w^2)^(1/2)+2*C(40)*D*Omega);
                (C(10)*S*rho*(u^2+v^2+w^2))/2,0,(C(11)*S*rho*(u^2+v^2+w^2))/2,0;
                0,(C(17)*S*rho*(u^2+v^2+w^2))/2,0,0;
                (C(23)*S*b*rho*(u^2+v^2+w^2))/2,0,(C(24)*S*b*rho*(u^2+v^2+w^2))/2,0;
                0,(C(30)*S*c*rho*(u^2+v^2+w^2))/2,0,0;
                (C(36)*S*b*rho*(u^2+v^2+w^2))/2,0,(C(37)*S*b*rho*(u^2+v^2+w^2))/2,0];
            
            dfdp = [(S*rho*(u^2+v^2+w^2))/2,(S*rho*atan2(w,u)*(u^2+v^2+w^2))/2,(S*c*q*rho*(u^2+v^2+w^2)^(1/2))/4,(S*de*rho*(u^2+v^2+w^2))/2,(S*rho*atan2(w,u)^2*(u^2+v^2+w^2))/2,(S*alphadot*rho*(u^2+v^2+w^2))/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,D^4*Omega^2*numprops*propefficiency*rho,D^3*Omega*numprops*propefficiency*rho*(u^2+v^2+w^2)^(1/2),D^2*numprops*propefficiency*rho*(u^2+v^2+w^2);
                0,0,0,0,0,0,(S*rho*asin(v/(u^2+v^2+w^2)^(1/2))*(u^2+v^2+w^2))/2,(S*b*p*rho*(u^2+v^2+w^2)^(1/2))/4,(S*b*r*rho*(u^2+v^2+w^2)^(1/2))/4,(S*da*rho*(u^2+v^2+w^2))/2,(S*dr*rho*(u^2+v^2+w^2))/2,(S*rho*asin(v/(u^2+v^2+w^2)^(1/2))^3*(u^2+v^2+w^2))/2,(S*betadot*rho*(u^2+v^2+w^2))/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                0,0,0,0,0,0,0,0,0,0,0,0,0,(S*rho*(u^2+v^2+w^2))/2,(S*rho*atan2(w,u)*(u^2+v^2+w^2))/2,(S*c*q*rho*(u^2+v^2+w^2)^(1/2))/4,(S*de*rho*(u^2+v^2+w^2))/2,(S*rho*atan2(w,u)^2*(u^2+v^2+w^2))/2,(S*alphadot*rho*(u^2+v^2+w^2))/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(S*b*rho*asin(v/(u^2+v^2+w^2)^(1/2))*(u^2+v^2+w^2))/2,(S*b^2*p*rho*(u^2+v^2+w^2)^(1/2))/4,(S*b^2*r*rho*(u^2+v^2+w^2)^(1/2))/4,(S*b*da*rho*(u^2+v^2+w^2))/2,(S*b*dr*rho*(u^2+v^2+w^2))/2,(S*b*rho*asin(v/(u^2+v^2+w^2)^(1/2))^3*(u^2+v^2+w^2))/2,(S*b*betadot*rho*(u^2+v^2+w^2))/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(S*c*rho*(u^2+v^2+w^2))/2,(S*c*rho*atan2(w,u)*(u^2+v^2+w^2))/2,(S*c^2*q*rho*(u^2+v^2+w^2)^(1/2))/4,(S*c*de*rho*(u^2+v^2+w^2))/2,(S*c*rho*atan2(w,u)^2*(u^2+v^2+w^2))/2,(S*alphadot*c*rho*(u^2+v^2+w^2))/2,0,0,0,0,0,0,0,0,0,0;
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,(S*b*rho*asin(v/(u^2+v^2+w^2)^(1/2))*(u^2+v^2+w^2))/2,(S*b^2*p*rho*(u^2+v^2+w^2)^(1/2))/4,(S*b^2*r*rho*(u^2+v^2+w^2)^(1/2))/4,(S*b*da*rho*(u^2+v^2+w^2))/2,(S*b*dr*rho*(u^2+v^2+w^2))/2,(S*b*rho*asin(v/(u^2+v^2+w^2)^(1/2))^3*(u^2+v^2+w^2))/2,(S*b*betadot*rho*(u^2+v^2+w^2))/2,0,0,0];

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
                newObj = FixedWing_NL_US(s.Aircraft,s.Name); 
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