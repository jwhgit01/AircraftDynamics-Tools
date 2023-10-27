classdef MultiRotor_BEMT_NASA < AerodynamicModel
%MultiRotor_BEMT_NASA

    properties
        % Nothing additional
    end

    methods
    
        function obj = MultiRotor_BEMT_NASA(Aircraft,Name)
            %MultiRotor_BEMT Construct an instance of this class
            %   Call the AerodynamicModel contstructor.
            obj@AerodynamicModel(Aircraft,Name);

            % Set the number of dynamic states and parameters
            obj.NumberOfStates = 0;
            obj.NumberOfParameters = 46;

            % % check to see if the CSV already exists
            % AircraftName = obj.Aircraft.Name;
            % AircraftPath = obj.Aircraft.WorkspaceDir;
            % if 0 ~= exist([AircraftPath filesep AircraftName '/AircraftData/' AircraftName '_Parameters_MultiRotor_BEMT.csv'],'file')
            %     warning(['A FixedWing_NL model already exists for ' AircraftName]);
            %     return
            % end
            % 
            % % Copy the template CSV
            % temp = which('AircraftDynamcis-Tools/lib/Parameters_MultiRotor_BEMT.csv');
            % [libpath,~,~] = fileparts(temp);
            % copyfile([libpath filesep 'Parameters_MultiRotor_BEMT.csv'],[AircraftPath filesep AircraftName '/AircraftData/' AircraftName '_Parameters_MultiRotor_BEMT.csv']);

        end % MultiRotor_BEMT

        function [Force,Moment] = Model(obj,~,vb,omega,Omega,Omegadot,Parameters,Constants)
            %Model

            % Parameters
            m = obj.Aircraft.Mass;
            h = obj.Aircraft.HubHeight;
            b = obj.Aircraft.ArmLength;
            Nr = obj.Aircraft.NumberOfRotors;
            R = obj.Aircraft.RotorRadius;
            M = obj.Aircraft.MixingMatrix;
            J = obj.Aircraft.MotorMomentOfInertia;

            % Constants
            rho = Constants.rho;
            g = Constants.g;
            rhopiR2 = rho*pi*R^2;
            nu0 = sqrt(m*g/(2*Nr*rhopiR2));
            
            % Airspeed and body velocity components 
            u = vb(1,1);
            v = vb(2,1);
            w = vb(3,1);
            V = sqrt(u^2 + v^2 + w^2);
            Vh = sqrt(u^2 + v^2);

            % Dynamic pressure
            qbar = 0.5*rho*V^2;

            % Aerodynamic angles
            if abs(u) > 1e-12
                % alpha = atan2(w,u);
                alpha = asin(w/V); % NASA Sim
            else
                alpha = 0;
            end
            if abs(V) > 1e-12
                % beta = asin(v/V);
                beta = atan2(v,u); % NASA Sim
            else
                beta = 0;
            end

            % Dimensional & nondimensional angular rates 
            p = omega(1,1);
            q = omega(2,1);
            r = omega(3,1);

            % Virtual actuator states
            delta = M*(Omega.^2);
            dt = delta(1,1);
            da = delta(2,1);
            de = delta(3,1);
            dr = delta(4,1);
            sqrtdelta = M*Omega;
            sqrtdt = sqrtdelta(1,1);
            sqrtda = sqrtdelta(2,1);
            sqrtde = sqrtdelta(3,1);
            sqrtdr = sqrtdelta(4,1);

            % Rotor acceleration 
            alpha_r = M(4,:)*Omegadot;

            % Parameters
            if isobject(Parameters)
                C = Parameters.Mean;
            else
                C = Parameters;
            end
            C_x_H_mu_b = C(1);
            C_x_H_mu_b_mu = C(2);
            C_x_H_mu_b_mu_w = C(3);
            C_x_H_mu_b_mu_0 = C(4);
            C_x_alphabeta = C(5);
            C_y_H_mu_b = C(6);
            C_y_H_mu_b_mu = C(7);
            C_y_H_mu_b_mu_w = C(8);
            C_y_H_mu_b_mu_0 = C(9);
            C_y_alphabeta = C(10);
            C_z_T_0 = C(11);
            C_z_T_mu = C(12);
            C_z_T_mu_w = C(13);
            C_z_T_mu_0 = C(14);
            C_z_T_mu2 = C(15);
            C_z_alpha = C(16);
            C_z_alphabeta = C(17);
            C_l_T_0 = C(18);
            C_l_T_mu = C(19);
            C_l_T_mu_w = C(20);
            C_l_T_mu_0 = C(21);
            C_l_H_mu_b = C(22);
            C_l_H_mu_b_mu = C(23);
            C_l_H_mu_b_mu_w = C(24);
            C_l_H_mu_b_mu_0 = C(25);
            C_l_R_mu_b = C(26);
            C_l_T_p = C(27);
            C_l_alphabeta = C(28);
            C_m_T_0 = C(29);
            C_m_T_mu = C(30);
            C_m_T_mu_w = C(31);
            C_m_T_mu_0 = C(32);
            C_m_H_mu_b = C(33);
            C_m_H_mu_b_mu = C(34);
            C_m_H_mu_b_mu_w = C(35);
            C_m_H_mu_b_mu_0 = C(36);
            C_m_R_mu_b = C(37);
            C_m_T_q = C(38);
            C_m_alphabeta = C(39);
            C_n_Q_0 = C(40);
            C_n_Q_mu = C(41);
            C_n_Q_mu_w = C(42);
            C_n_Q_mu_0 = C(43);
            C_n_H_mu_b = C(44);
            C_n_T_r = C(45);
            C_n_alphabeta = C(46);

            % Forces and moment components
            %
            X = rhopiR2*u*(-C_x_H_mu_b*R*sqrtdt...
                           -C_x_H_mu_b_mu*Nr*Vh...
                           +C_x_H_mu_b_mu_w*Nr*w...
                           -C_x_H_mu_b_mu_0*Nr*nu0)...
                +qbar*C_x_alphabeta*cos(alpha)*cos(beta);
            %
            Y = rhopiR2*v*(-C_y_H_mu_b*R*sqrtdt...
                           -C_y_H_mu_b_mu*Nr*Vh...
                           +C_y_H_mu_b_mu_w*Nr*w...
                           -C_y_H_mu_b_mu_0*Nr*nu0)...
                +qbar*C_y_alphabeta*cos(alpha)*sin(beta);
            %
            Z = rhopiR2*(-C_z_T_0*R^2*dt...
                         +C_z_T_mu*R*Vh*sqrtdt...
                         +C_z_T_mu_w*R*w*sqrtdt...
                         +C_z_T_mu_0*R*nu0*sqrtdt...
                         -C_z_T_mu2*Nr*Vh^2)...
                +qbar*(C_z_alpha*sin(alpha)...
                          +C_z_alphabeta*sin(alpha)*cos(Nr*beta));
            %
            L = q*J*sqrtdr + rhopiR2*(+C_l_T_0*R^2*da...
                         -C_l_T_mu*R*Vh*sqrtda...
                         +C_l_T_mu_w*R*w*sqrtda...
                         -C_l_T_mu_0*R*nu0*sqrtda...
                         -C_l_H_mu_b*h*R*v*sqrtdt...
                         -C_l_H_mu_b_mu*Nr*h*v*Vh...
                         +C_l_H_mu_b_mu_w*Nr*h*v*w...
                         -C_l_H_mu_b_mu_0*Nr*h*v*nu0...
                         +C_l_R_mu_b*R^2*u*sqrtdr...
                         +C_l_T_p*R*b*sqrtdt*p)...
                +qbar*C_l_alphabeta*cos(alpha)*sin(beta);
            %
            M = -p*J*sqrtdr + rhopiR2*(+C_m_T_0*R^2*de...
                         -C_m_T_mu*R*Vh*sqrtde...
                         +C_m_T_mu_w*R*w*sqrtde...
                         -C_m_T_mu_0*R*nu0*sqrtde...
                         +C_m_H_mu_b*h*R*u*sqrtdt...
                         +C_m_H_mu_b_mu*Nr*h*u*Vh...
                         -C_m_H_mu_b_mu_w*Nr*h*u*w...
                         +C_m_H_mu_b_mu_0*Nr*h*u*nu0...
                         +C_m_R_mu_b*R^2*v*sqrtdr...
                         +C_m_T_q*R*b*sqrtdt*q)...
                +qbar*C_m_alphabeta*sin(2*alpha)*cos(beta);
            %
            N = J*alpha_r + rhopiR2*(-C_n_Q_0*R^3*dr...
                         -C_n_Q_mu*R^2*Vh*sqrtdr...
                         +C_n_Q_mu_w*R^2*w*sqrtdr...
                         -C_n_Q_mu_0*R^2*nu0*sqrtdr...
                         -C_n_H_mu_b*R*(u*sqrtda+v*sqrtde)...
                         +C_n_T_r*R*b*sqrtdt*r)...
                +qbar*C_n_alphabeta*cos(alpha)*sin(2*beta);

            % Force and moment vectors
            Force = [X;Y;Z];
            Moment = [L;M;N];
            
        end % Model

        function detadt = AugmentedDynamics(obj,~,~,~,~,~,~)
            detadt = zeros(0,1);
        end % AugmentedDynamics

        function [dfdv,dfdomega,dfdu,dfdp,dfdeta] = Jacobian(obj,~,vb,omega,Omega,Omegadot,Parameters,Constants)
            %Jacobian

            % Parameters
            m = obj.Aircraft.Mass;
            h = obj.Aircraft.HubHeight;
            b = obj.Aircraft.ArmLength;
            Nr = obj.Aircraft.NumberOfRotors;
            R = obj.Aircraft.RotorRadius;
            M = obj.Aircraft.MixingMatrix;
            J = obj.Aircraft.MotorMomentOfInertia;

            % Constants
            rho = Constants.rho;
            g = Constants.g;
            rhopiR2 = rho*pi*R^2;
            nu0 = sqrt(m*g/(2*Nr*rhopiR2));
            
            % Airspeed and body velocity components 
            u = vb(1,1);
            v = vb(2,1);
            w = vb(3,1);
            V = sqrt(u^2 + v^2 + w^2);
            Vh = sqrt(u^2 + v^2);

            % Dynamic pressure
            qbar = 0.5*rho*V^2;

            % Aerodynamic angles
            if abs(u) > 1e-12
                % alpha = atan2(w,u);
                alpha = asin(w/V); % NASA Sim
            else
                alpha = 0;
            end
            if abs(V) > 1e-12
                % beta = asin(v/V);
                beta = atan2(v,u); % NASA Sim
            else
                beta = 0;
            end

            % Angular rates 
            p = omega(1,1);
            q = omega(2,1);
            r = omega(3,1);

            % Virtual actuator states
            delta = M*(Omega.^2);
            dt = delta(1,1);
            da = delta(2,1);
            de = delta(3,1);
            dr = delta(4,1);
            sqrtdelta = M*Omega;
            sqrtdt = sqrtdelta(1,1);
            sqrtda = sqrtdelta(2,1);
            sqrtde = sqrtdelta(3,1);
            sqrtdr = sqrtdelta(4,1);

            % Rotor acceleration 
            % alpha_r = M(4,:)*Omegadot;

            % Parameters
            if isobject(Parameters)
                C = Parameters.Mean;
            else
                C = Parameters;
            end

            % The following were computed using symbolic computation on 9/15/2023.
            dfdv = [(C(5)*rho*u^2*(1-v^2/(u^2+v^2+w^2))^(1/2))/(u^2+w^2)^(1/2)-(C(2)*Nr*rhopiR2*u^2)/(u^2+v^2)^(1/2)-rhopiR2*(C(2)*Nr*(u^2+v^2)^(1/2)+C(4)*Nr*nu0+C(1)*R*sqrtdt-C(3)*Nr*w)+(C(5)*rho*(1-v^2/(u^2+v^2+w^2))^(1/2)*(u^2+v^2+w^2))/(2*(u^2+w^2)^(1/2))-(C(5)*rho*u^2*(1-v^2/(u^2+v^2+w^2))^(1/2)*(u^2+v^2+w^2))/(2*(u^2+w^2)^(3/2))+(C(5)*rho*u^2*v^2)/(2*(u^2+w^2)^(1/2)*(1-v^2/(u^2+v^2+w^2))^(1/2)*(u^2+v^2+w^2)),(C(5)*rho*u*v*(1-v^2/(u^2+v^2+w^2))^(1/2))/(u^2+w^2)^(1/2)-(C(2)*Nr*rhopiR2*u*v)/(u^2+v^2)^(1/2)-(C(5)*rho*u*((2*v)/(u^2+v^2+w^2)-(2*v^3)/(u^2+v^2+w^2)^2)*(u^2+v^2+w^2))/(4*(u^2+w^2)^(1/2)*(1-v^2/(u^2+v^2+w^2))^(1/2)),C(3)*Nr*rhopiR2*u+(C(5)*rho*u*w*(1-v^2/(u^2+v^2+w^2))^(1/2))/(u^2+w^2)^(1/2)-(C(5)*rho*u*w*(1-v^2/(u^2+v^2+w^2))^(1/2)*(u^2+v^2+w^2))/(2*(u^2+w^2)^(3/2))+(C(5)*rho*u*v^2*w)/(2*(u^2+w^2)^(1/2)*(1-v^2/(u^2+v^2+w^2))^(1/2)*(u^2+v^2+w^2));
                (C(10)*rho*v*(u^2+v^2+w^2)^(1/2))/(2*(u^2+w^2)^(1/2))-(C(7)*Nr*rhopiR2*u*v)/(u^2+v^2)^(1/2)+(C(10)*rho*u^2*v)/(2*(u^2+w^2)^(1/2)*(u^2+v^2+w^2)^(1/2))-(C(10)*rho*u^2*v*(u^2+v^2+w^2)^(1/2))/(2*(u^2+w^2)^(3/2)),(C(10)*rho*u*(u^2+v^2+w^2)^(1/2))/(2*(u^2+w^2)^(1/2))-(C(7)*Nr*rhopiR2*v^2)/(u^2+v^2)^(1/2)-rhopiR2*(C(7)*Nr*(u^2+v^2)^(1/2)+C(9)*Nr*nu0+C(6)*R*sqrtdt-C(8)*Nr*w)+(C(10)*rho*u*v^2)/(2*(u^2+w^2)^(1/2)*(u^2+v^2+w^2)^(1/2)),C(8)*Nr*rhopiR2*v+(C(10)*rho*u*v*w)/(2*(u^2+w^2)^(1/2)*(u^2+v^2+w^2)^(1/2))-(C(10)*rho*u*v*w*(u^2+v^2+w^2)^(1/2))/(2*(u^2+w^2)^(3/2));
                rho*u*((C(16)*w)/(u^2+w^2)^(1/2)+(C(17)*w*cos(Nr*asin(v/(u^2+v^2+w^2)^(1/2))))/(u^2+w^2)^(1/2))-(rho*((C(16)*u*w)/(u^2+w^2)^(3/2)+(C(17)*u*w*cos(Nr*asin(v/(u^2+v^2+w^2)^(1/2))))/(u^2+w^2)^(3/2)-(C(17)*Nr*u*v*w*sin(Nr*asin(v/(u^2+v^2+w^2)^(1/2))))/((u^2+w^2)^(1/2)*(1-v^2/(u^2+v^2+w^2))^(1/2)*(u^2+v^2+w^2)^(3/2)))*(u^2+v^2+w^2))/2-rhopiR2*(2*C(15)*Nr*u-(C(12)*R*sqrtdt*u)/(u^2+v^2)^(1/2)),rho*v*((C(16)*w)/(u^2+w^2)^(1/2)+(C(17)*w*cos(Nr*asin(v/(u^2+v^2+w^2)^(1/2))))/(u^2+w^2)^(1/2))-rhopiR2*(2*C(15)*Nr*v-(C(12)*R*sqrtdt*v)/(u^2+v^2)^(1/2))-(C(17)*Nr*rho*w*sin(Nr*asin(v/(u^2+v^2+w^2)^(1/2)))*(1/(u^2+v^2+w^2)^(1/2)-v^2/(u^2+v^2+w^2)^(3/2))*(u^2+v^2+w^2))/(2*(u^2+w^2)^(1/2)*(1-v^2/(u^2+v^2+w^2))^(1/2)),C(13)*R*rhopiR2*sqrtdt+rho*w*((C(16)*w)/(u^2+w^2)^(1/2)+(C(17)*w*cos(Nr*asin(v/(u^2+v^2+w^2)^(1/2))))/(u^2+w^2)^(1/2))+(rho*(u^2+v^2+w^2)*(C(16)/(u^2+w^2)^(1/2)+(C(17)*cos(Nr*asin(v/(u^2+v^2+w^2)^(1/2))))/(u^2+w^2)^(1/2)-(C(16)*w^2)/(u^2+w^2)^(3/2)-(C(17)*w^2*cos(Nr*asin(v/(u^2+v^2+w^2)^(1/2))))/(u^2+w^2)^(3/2)+(C(17)*Nr*v*w^2*sin(Nr*asin(v/(u^2+v^2+w^2)^(1/2))))/((u^2+w^2)^(1/2)*(1-v^2/(u^2+v^2+w^2))^(1/2)*(u^2+v^2+w^2)^(3/2))))/2;
                (C(28)*rho*v*(u^2+v^2+w^2)^(1/2))/(2*(u^2+w^2)^(1/2))-rhopiR2*((C(19)*R*sqrtda*u)/(u^2+v^2)^(1/2)-C(26)*R^2*sqrtdr+(C(23)*Nr*h*u*v)/(u^2+v^2)^(1/2))+(C(28)*rho*u^2*v)/(2*(u^2+w^2)^(1/2)*(u^2+v^2+w^2)^(1/2))-(C(28)*rho*u^2*v*(u^2+v^2+w^2)^(1/2))/(2*(u^2+w^2)^(3/2)),(C(28)*rho*u*(u^2+v^2+w^2)^(1/2))/(2*(u^2+w^2)^(1/2))-rhopiR2*(C(23)*Nr*h*(u^2+v^2)^(1/2)+C(25)*Nr*h*nu0+C(22)*R*h*sqrtdt-C(24)*Nr*h*w+(C(19)*R*sqrtda*v)/(u^2+v^2)^(1/2)+(C(23)*Nr*h*v^2)/(u^2+v^2)^(1/2))+(C(28)*rho*u*v^2)/(2*(u^2+w^2)^(1/2)*(u^2+v^2+w^2)^(1/2)),rhopiR2*(C(20)*R*sqrtda+C(24)*Nr*h*v)+(C(28)*rho*u*v*w)/(2*(u^2+w^2)^(1/2)*(u^2+v^2+w^2)^(1/2))-(C(28)*rho*u*v*w*(u^2+v^2+w^2)^(1/2))/(2*(u^2+w^2)^(3/2));
                rhopiR2*(C(34)*Nr*h*(u^2+v^2)^(1/2)+C(36)*Nr*h*nu0+C(33)*R*h*sqrtdt-C(35)*Nr*h*w-(C(30)*R*sqrtde*u)/(u^2+v^2)^(1/2)+(C(34)*Nr*h*u^2)/(u^2+v^2)^(1/2))+C(39)*rho*u*sin(2*atan2(w,u))*(1-v^2/(u^2+v^2+w^2))^(1/2)-(C(39)*rho*w*cos(2*atan2(w,u))*(1-v^2/(u^2+v^2+w^2))^(1/2)*(u^2+v^2+w^2))/(u^2+w^2)+(C(39)*rho*u*v^2*sin(2*atan2(w,u)))/(2*(1-v^2/(u^2+v^2+w^2))^(1/2)*(u^2+v^2+w^2)),rhopiR2*(C(37)*R^2*sqrtdr-(C(30)*R*sqrtde*v)/(u^2+v^2)^(1/2)+(C(34)*Nr*h*u*v)/(u^2+v^2)^(1/2))+C(39)*rho*v*sin(2*atan2(w,u))*(1-v^2/(u^2+v^2+w^2))^(1/2)-(C(39)*rho*sin(2*atan2(w,u))*((2*v)/(u^2+v^2+w^2)-(2*v^3)/(u^2+v^2+w^2)^2)*(u^2+v^2+w^2))/(4*(1-v^2/(u^2+v^2+w^2))^(1/2)),rhopiR2*(C(31)*R*sqrtde-C(35)*Nr*h*u)+C(39)*rho*w*sin(2*atan2(w,u))*(1-v^2/(u^2+v^2+w^2))^(1/2)+(C(39)*rho*v^2*w*sin(2*atan2(w,u)))/(2*(1-v^2/(u^2+v^2+w^2))^(1/2)*(u^2+v^2+w^2))+(C(39)*rho*u*cos(2*atan2(w,u))*(1-v^2/(u^2+v^2+w^2))^(1/2)*(u^2+v^2+w^2))/(u^2+w^2);
                (C(46)*rho*u^2*sin(2*asin(v/(u^2+v^2+w^2)^(1/2))))/(u^2+w^2)^(1/2)-rhopiR2*(C(44)*R*sqrtda+(C(41)*R^2*sqrtdr*u)/(u^2+v^2)^(1/2))+(C(46)*rho*sin(2*asin(v/(u^2+v^2+w^2)^(1/2)))*(u^2+v^2+w^2))/(2*(u^2+w^2)^(1/2))-(C(46)*rho*u^2*sin(2*asin(v/(u^2+v^2+w^2)^(1/2)))*(u^2+v^2+w^2))/(2*(u^2+w^2)^(3/2))-(C(46)*rho*u^2*v*cos(2*asin(v/(u^2+v^2+w^2)^(1/2))))/((u^2+w^2)^(1/2)*(1-v^2/(u^2+v^2+w^2))^(1/2)*(u^2+v^2+w^2)^(1/2)),(C(46)*rho*u*v*sin(2*asin(v/(u^2+v^2+w^2)^(1/2))))/(u^2+w^2)^(1/2)-rhopiR2*(C(44)*R*sqrtde+(C(41)*R^2*sqrtdr*v)/(u^2+v^2)^(1/2))+(C(46)*rho*u*cos(2*asin(v/(u^2+v^2+w^2)^(1/2)))*(1/(u^2+v^2+w^2)^(1/2)-v^2/(u^2+v^2+w^2)^(3/2))*(u^2+v^2+w^2))/((u^2+w^2)^(1/2)*(1-v^2/(u^2+v^2+w^2))^(1/2)),C(42)*R^2*rhopiR2*sqrtdr+(C(46)*rho*u*w*sin(2*asin(v/(u^2+v^2+w^2)^(1/2))))/(u^2+w^2)^(1/2)-(C(46)*rho*u*w*sin(2*asin(v/(u^2+v^2+w^2)^(1/2)))*(u^2+v^2+w^2))/(2*(u^2+w^2)^(3/2))-(C(46)*rho*u*v*w*cos(2*asin(v/(u^2+v^2+w^2)^(1/2))))/((u^2+w^2)^(1/2)*(1-v^2/(u^2+v^2+w^2))^(1/2)*(u^2+v^2+w^2)^(1/2))];
            %
            dfdomega = [0,0,0;
                0,0,0;
                0,0,0;
                C(27)*R*b*rhopiR2*sqrtdt,J*sqrtdr,0;
                -J*sqrtdr,C(38)*R*b*rhopiR2*sqrtdt,0;
                0,0,C(45)*R*b*rhopiR2*sqrtdt];
            %
            dfdu = [0,0,0,0,-C(1)*R*rhopiR2*u,0,0,0,0;
                0,0,0,0,-C(6)*R*rhopiR2*v,0,0,0,0;
                -C(11)*R^2*rhopiR2,0,0,0,rhopiR2*(C(12)*R*Vh+C(14)*R*nu0+C(13)*R*w),0,0,0,0;
                0,C(18)*R^2*rhopiR2,0,0,rhopiR2*(C(27)*R*b*p-C(22)*R*h*v),-rhopiR2*(C(19)*R*Vh+C(21)*R*nu0-C(20)*R*w),0,C(26)*rhopiR2*u*R^2+J*q,0;
                0,0,C(29)*R^2*rhopiR2,0,rhopiR2*(C(38)*R*b*q+C(33)*R*h*u),0,-rhopiR2*(C(30)*R*Vh+C(32)*R*nu0-C(31)*R*w),C(37)*rhopiR2*v*R^2-J*p,0;
                0,0,0,-C(40)*R^3*rhopiR2,C(45)*R*b*r*rhopiR2,-C(44)*R*rhopiR2*u,-C(44)*R*rhopiR2*v,-rhopiR2*(C(41)*R^2*Vh+C(43)*R^2*nu0-C(42)*R^2*w),J];
            %
            dfdp = [-R*rhopiR2*sqrtdt*u,-Nr*Vh*rhopiR2*u,Nr*rhopiR2*u*w,-Nr*nu0*rhopiR2*u,qbar*cos(alpha)*cos(beta),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                0,0,0,0,0,-R*rhopiR2*sqrtdt*v,-Nr*Vh*rhopiR2*v,Nr*rhopiR2*v*w,-Nr*nu0*rhopiR2*v,qbar*cos(alpha)*sin(beta),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                0,0,0,0,0,0,0,0,0,0,-R^2*dt*rhopiR2,R*Vh*rhopiR2*sqrtdt,R*rhopiR2*sqrtdt*w,R*nu0*rhopiR2*sqrtdt,-Nr*Vh^2*rhopiR2,qbar*sin(alpha),qbar*cos(Nr*beta)*sin(alpha),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,R^2*da*rhopiR2,-R*Vh*rhopiR2*sqrtda,R*rhopiR2*sqrtda*w,-R*nu0*rhopiR2*sqrtda,-R*h*rhopiR2*sqrtdt*v,-Nr*Vh*h*rhopiR2*v,Nr*h*rhopiR2*v*w,-Nr*h*nu0*rhopiR2*v,R^2*rhopiR2*sqrtdr*u,R*b*p*rhopiR2*sqrtdt,qbar*cos(alpha)*sin(beta),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,R^2*de*rhopiR2,-R*Vh*rhopiR2*sqrtde,R*rhopiR2*sqrtde*w,-R*nu0*rhopiR2*sqrtde,R*h*rhopiR2*sqrtdt*u,Nr*Vh*h*rhopiR2*u,-Nr*h*rhopiR2*u*w,Nr*h*nu0*rhopiR2*u,R^2*rhopiR2*sqrtdr*v,R*b*q*rhopiR2*sqrtdt,qbar*sin(2*alpha)*cos(beta),0,0,0,0,0,0,0;
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-R^3*dr*rhopiR2,-R^2*Vh*rhopiR2*sqrtdr,R^2*rhopiR2*sqrtdr*w,-R^2*nu0*rhopiR2*sqrtdr,-R*rhopiR2*(sqrtda*u+sqrtde*v),R*b*r*rhopiR2*sqrtdt,qbar*sin(2*beta)*cos(alpha)];
            %
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
                newObj = MultiRotor_BEMT_NASA(s.Aircraft,s.Name); 
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