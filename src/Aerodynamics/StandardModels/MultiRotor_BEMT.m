classdef MultiRotor_BEMT < AerodynamicModel
%MultiRotor_BEMT

    properties
        % Nothing additional
    end

    methods
    
        function obj = MultiRotor_BEMT(Aircraft,Name)
            %MultiRotor_BEMT Construct an instance of this class
            %   Call the AerodynamicModel contstructor.
            obj@AerodynamicModel(Aircraft,Name);

            % Set the number of dynamic states and parameters
            obj.NumberOfStates = 0;
            obj.NumberOfParameters = 44;

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
            Vh = sqrt(u^2 + v^2);

            % Dimensional & nondimensional angular rates 
            p = omega(1,1);
            q = omega(2,1);
            r = omega(3,1);
            phat = p*b/nu0;
            qhat = q*b/nu0;
            rhat = r*b/nu0;

            % Virtual actuator states
            delta2 = M*(Omega.^2);
            d2t = delta2(1,1);
            d2a = delta2(2,1);
            d2e = delta2(3,1);
            d2r = delta2(4,1);
            delta = M*Omega;
            dt = delta(1,1);
            da = delta(2,1);
            de = delta(3,1);
            dr = delta(4,1);

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
            C_x_u = C(5);
            C_y_H_mu_b = C(6);
            C_y_H_mu_b_mu = C(7);
            C_y_H_mu_b_mu_w = C(8);
            C_y_H_mu_b_mu_0 = C(9);
            C_y_v = C(10);
            C_z_T_0 = C(11);
            C_z_T_mu = C(12);
            C_z_T_mu_w = C(13);
            C_z_T_mu_0 = C(14);
            C_z_T_mu2 = C(15);
            C_z_w = C(16);
            C_l_T_0 = C(17);
            C_l_T_mu = C(18);
            C_l_T_mu_w = C(19);
            C_l_T_mu_0 = C(20);
            C_l_H_mu_b = C(21);
            C_l_H_mu_b_mu = C(22);
            C_l_H_mu_b_mu_w = C(23);
            C_l_H_mu_b_mu_0 = C(24);
            C_l_R_mu_b = C(25);
            C_l_p = C(26);
            C_m_T_0 = C(27);
            C_m_T_mu = C(28);
            C_m_T_mu_w = C(29);
            C_m_T_mu_0 = C(30);
            C_m_H_mu_b = C(31);
            C_m_H_mu_b_mu = C(32);
            C_m_H_mu_b_mu_w = C(33);
            C_m_H_mu_b_mu_0 = C(34);
            C_m_R_mu_b = C(35);
            C_m_q = C(36);
            C_n_Q_0 = C(37);
            C_n_Q_mu = C(38);
            C_n_Q_mu_w = C(39);
            C_n_Q_mu_0 = C(40);
            C_n_Q_mu_w2 = C(41);
            C_n_H_mu_b = C(42);
            C_n_H_mu_b_mu_w = C(43);
            C_n_r = C(44);

            % Forces and moment components
            %
            X = rhopiR2*Nr*u*(-C_x_H_mu_b*R*dt...
                           -C_x_H_mu_b_mu*Vh...
                           +C_x_H_mu_b_mu_w*w...
                           -C_x_H_mu_b_mu_0*nu0)...
                +0.5*rho*u^2*b*h*C_x_u*sign(u);
            %
            Y = rhopiR2*Nr*v*(-C_y_H_mu_b*R*dt...
                           -C_y_H_mu_b_mu*Vh...
                           +C_y_H_mu_b_mu_w*w...
                           -C_y_H_mu_b_mu_0*nu0)...
                +0.5*rho*v^2*b*h*C_y_v*sign(v);
            %
            Z = rhopiR2*(-C_z_T_0*Nr*R^2*d2t...
                         +C_z_T_mu*Nr*R*Vh*dt...
                         +C_z_T_mu_w*R*(p*da + q*de - Nr*w*dt)...
                         +C_z_T_mu_0*Nr*R*nu0*dt...
                         -C_z_T_mu2*Nr*Vh^2)...
                +0.5*rho*w^2*b^2*C_z_w*sign(w);
            %
            L = q*J*dr + rhopiR2*(+C_l_T_0*R^2*d2a...
                         -C_l_T_mu*R*Vh*da...
                         +C_l_T_mu_w*R*(w*da - Nr*b^2*p*dt + b^2*q*dr)...
                         -C_l_T_mu_0*R*nu0*da...
                         -C_l_H_mu_b*Nr*R*h*v*dt...
                         -C_l_H_mu_b_mu*Nr*h*v*Vh...
                         +C_l_H_mu_b_mu_w*Nr*h*v*w...
                         -C_l_H_mu_b_mu_0*Nr*h*v*nu0...
                         +C_l_R_mu_b*R^2*u*dr)...
                +0.5*rho*(v^2+w^2)*b^2*h*C_l_p*phat;
            %
            M = -p*J*dr + rhopiR2*(+C_m_T_0*R^2*d2e...
                         -C_m_T_mu*R*Vh*de...
                         +C_m_T_mu_w*R*(w*de - Nr*b^2*q*dt + b^2*p*dr)...
                         -C_m_T_mu_0*R*nu0*de...
                         +C_m_H_mu_b*Nr*R*h*u*dt...
                         +C_m_H_mu_b_mu*Nr*h*u*Vh...
                         -C_m_H_mu_b_mu_w*Nr*h*u*w...
                         +C_m_H_mu_b_mu_0*Nr*h*u*nu0...
                         +C_m_R_mu_b*R^2*v*dr)...
                +0.5*rho*(u^2+w^2)*b^2*h*C_m_q*qhat;
            %
            N = J*alpha_r + rhopiR2*(-C_n_Q_0*R^3*d2r...
                         -C_n_Q_mu*R^2*Vh*dr...
                         +C_n_Q_mu_w*R^2*(w*dr + p*de + q*da)...
                         -C_n_Q_mu_0*R^2*nu0*dr...
                         -C_n_Q_mu_w2*2*Nr*R*b^2*p*q...
                         -C_n_H_mu_b*R*(u*da+v*de)...
                         -C_n_H_mu_b_mu_w*Nr*b^2*(u*p + v*q))...
                +0.5*rho*(u^2+v^2)*b^2*h*C_n_r*rhat;

            % Force and moment vectors
            Force = [X;Y;Z];
            Moment = [L;M;N];
            
        end % Model

        function [H,X,info] = RegressorMatrix(obj,vb,omega,Omega,Constants)
            %RegressorMatrix
            %
            % This function computes the force and moment regressors for a blade
            % element momentum theory (BEMT) model of a multirotor aircraft. This model
            % is linearly parameterized and has the following measurement model for the
            % k'th sample:
            % 
            %               y(k) = Hk(vb(k),omega(k),Omega(k))*theta            (1)
            %
            % where y(k) = [F(k); M(k)-(motor intertial moments)], 
            % Hk: R^(6+Nr) --> R^{6 x np}, theta is the np x 1 parameter vector. Note
            % that this function does not treat the motor moment of inertia as an 
            % unknown parameter. Instead its associated moments are considered part of
            % the model output, y.
            %
            % If model parameters are assumed to be shared among axes, the regressors
            % in H correspond to the parameter vector
            %
            % theta = [C_H_mub, C_H_mub_mu0, C_H_mub_mu, C_H_mub_muw, C_T_0,...
            %          C_T_mu0, C_T_mu, C_T_mu2, C_T_muw, C_R_mub, C_Q_0, C_Q_mu0,...
            %          C_Q_mu, C_Q_muw, (C_Q_muw2), C_x_u, C_y_v, C_z_w, C_l_p,...
            %          C_m_q, C_n_r]'
            %
            % where parameters in parentheses may vanish in certain configurations. In
            % this case, the N-sample time history of regressors is given as
            %
            %                      H = [H(1);H(2);...;H(N)]                     (2)
            %
            % If model parameters are assumed to not be shared among axes, we consider
            % the following measurment model instead. For each axis, the measurement
            % model for the k'th sample is
            %
            %               y(k) = x'(vb(k),omega(k),Omega(k))*theta            (3)
            %
            % where x is a np x 1 regressor vector. For this model, the regressors
            % in x correspond to each of the following parameter vectors
            %
            %   Fx: [C_X_H_mub, C_X_H_mub_mu0, C_X_H_mub_mu, C_X_H_mub_muw, C_x_u]'
            %   Fy: [C_Y_H_mub, C_Y_H_mub_mu0, C_Y_H_mub_mu, C_Y_H_mub_muw, C_y_v]'
            %   Fz: [C_Z_T_0, C_Z_T_mu0, C_Z_T_mu, C_Z_T_mu2, C_Z_T_muw, C_z_w]'
            %   Mx: [C_l_R_mub, C_l_T_0, C_l_T_mu0, C_l_T_mu, C_l_T_muw, C_l_H_mub, ...
            %        C_l_H_mub_mu0, C_l_H_mub_mu, C_l_H_mub_muw, C_l_p]'
            %   My: [C_m_R_mub, C_m_T_0, C_m_T_mu0, C_m_T_mu, C_m_T_muw, C_m_H_mub, ...
            %        C_m_H_mub_mu0, C_m_H_mub_mu, C_m_H_mub_muw, C_m_q]'
            %   Mz: [C_n_Q_0, C_n_Q_mu0, C_n_Q_mu, C_n_Q_muw, (C_n_Q_muw2), ...
            %        C_n_H_mub, C_n_H_mub_muw, C_n_r]'
            %
            % In this case, the N-sample time history of regressors for each axis is
            %
            %                      X = [x'(1);x'(2);...;x'(N)]                  (4)
            %
            % Inputs:
            %
            %   vb          A 3 x N array containing the relative velocity vector
            %               expressed in the body frame.
            %
            %   omega       A 3 x N array containing the relative angular velocity
            %               vector of the body frame.
            %
            %   Omega       A Nr x N array containing rotor speeds (non-negative).
            %
            %   Constants   A struct containing the following enviromental and aircraft
            %               constants:
            %                   g               gravitational acceleration
            %                   rho             air density
            %                   Mass            aircraft mass
            %                   HubHeight       height of rotor hub above aircraft CG
            %                   ArmLength       rotor arm length
            %                   NumberOfRotors  number of rotors
            %                   RotorRadius     propeller radius
            %                   MixingMatrix    aircraft model mixing matrix
            %                   PlusConfig      logical indicating whether rotors are
            %                                   in a "+" configuration (a rotor is
            %                                   located along the b1 axis)
            %                   MotorMoI        motor/propeller moment of inertia about
            %                                   z-axis
            %                   A12, A13, A23   reference areas in the b1-b2, b1-b3,
            %                                   and b2-b3 planes, respectively
            %  
            % Outputs:
            %
            %   H       The N*6 x np matrix defined in (2) that contains the model
            %           regressors for the measurment model in (1).
            %
            %   X       A 6 x 1 cell array with the i'th cell containing a N x npi
            %           regressor matrix defined in (4) for the i'th foce/moment axis
            %           according to the model in (3).
            %
            %   info    TBD
            %
            
            % Aircraft Constants
            m = Constants.Mass;
            h = Constants.HubHeight;
            l = Constants.ArmLength;
            Nr = Constants.NumberOfRotors;
            R = Constants.RotorRadius;
            M = Constants.MixingMatrix;
            PlusConfig = Constants.PlusConfig;
            % J = Constants.MotorMoI; % not currently used
            % A12 = Constants.A12; % not currently used
            A13 = Constants.A13;
            A23 = Constants.A23;
            
            % Enviromental Constants
            rho = Constants.rho;
            g = Constants.g;
            
            % Compute commonly used constants for efficiency
            R2 = R^2;
            R3 = R^3;
            rhopiR2 = rho*pi*R2;
            nu0 = sqrt(m*g/(2*Nr*rhopiR2));
            onehalfl2 = 0.5*l^2;
            info.nu0 = nu0;
            
            % Check the number of samples
            N = size(vb,2);
            if size(omega,2) ~= N || size(Omega,2) ~= N
                error('Explanatory variable dimensions do not agree.');
            end
            
            % Number of parameters per axis
            np.Fx = 5;
            np.Fy = 5;
            np.Fz = 6;
            np.Mx = 10;
            np.My = 10;
            if Nr == 4
                np.Mz = 8;
                np.Total = 21;
            else
                np.Mz = 7;
                np.Total = 20;
            end
            info.NumberOfParameters = np;
            
            % Initialize results
            H = zeros(6*N,np.Total);
            X = cell(6,1);
            for ii = 1:6
                X{ii} = zeros(N,1);
            end
            
            % Loop through the N samples to construct H and X
            for k = 1:N
            
                % Get the explanatory variables for the k'th sample
                u = vb(1,k);
                v = vb(2,k);
                w = vb(3,k);
                Vh = sqrt(u^2 + v^2);
                p = omega(1,k);
                q = omega(2,k);
                r = omega(3,k);
                delta2 = M*(Omega(:,k).^2);
                d2t = delta2(1,k);
                d2a = delta2(2,k);
                d2e = delta2(3,k);
                d2r = delta2(4,k);
                delta = M*Omega(:,k);
                dt = delta(1,k);
                da = delta(2,k);
                de = delta(3,k);
                dr = delta(4,k);
            
                % Configuration-dependent regressors
                if Nr == 4 && PlusConfig
                    Delta_omegadelta_w_l = onehalfl2*p*(Nr*dt - dr);
                    Delta_omegadelta_w_m = onehalfl2*q*(Nr*dt + dr);
                    Delta_omegadelta_w_n = p*da - q*de;
                    Delta_omega_w2_n = onehalfl2*Nr*(p^2 - q^2);
                elseif Nr == 4 && ~PlusConfig
                    Delta_omegadelta_w_l = onehalfl2*(Nr*p*dt + q*dr);
                    Delta_omegadelta_w_m = onehalfl2*(Nr*q*dt + p*dr);
                    Delta_omegadelta_w_n = -p*de - q*da;
                    Delta_omega_w2_n = -2*onehalfl2*Nr*p*q;
                elseif Nr >= 6
                    Delta_omegadelta_w_l = onehalfl2*p*Nr*dt;
                    Delta_omegadelta_w_m = onehalfl2*q*Nr*dt;
                    Delta_omegadelta_w_n = 0;
                    Delta_omega_w2_n = 0;
                else
                    error('Invalid configuration');
                end
                
                % Fx
                xk_X_H_mub = -rhopiR2*Nr*u*dt;
                xk_X_H_mub_mu0 = -rhopiR2*Nr*u*nu0;
                xk_X_H_mub_mu = -rhopiR2*Nr*u*Vh;
                xk_X_H_mub_muw = +rhopiR2*Nr*u*w;
                xk_x_u = +0.5*rho*u^2*A23*sign(u);
            
                % Fy
                xk_Y_H_mub = -rhopiR2*Nr*v*dt;
                xk_Y_H_mub_mu0 = -rhopiR2*Nr*v*nu0;
                xk_Y_H_mub_mu = -rhopiR2*Nr*v*Vh;
                xk_Y_H_mub_muw = +rhopiR2*Nr*v*w;
                xk_y_v = +0.5*rho*v^2*A13*sign(v);
                
                % Fz
                xk_Z_T_0 = -rhopiR2*R2*Nr*d2t;
                xk_Z_T_mu0 = +rhopiR2*R*Nr*nu0*dt;
                xk_Z_T_mu = +rhopiR2*R*Nr*Vh*dt;
                xk_Z_T_mu2 = -rhopiR2*Vh^2;
                xk_Z_T_muw = -rhopiR2*R*(Nr*w*dt - p*da - q*de);
                xk_z_w = +0.5*rho*w^2*A23*sign(w);
            
                % Mx
                xk_l_R_mub = +rhopiR2*R2*u*dr;
                xk_l_T_0 = +rhopiR2*R2*d2a;
                xk_l_T_mu0 = -rhopiR2*R*nu0*da;
                xk_l_T_mu = -rhopiR2*R*Vh*da;
                xk_l_T_muw = +rhopiR2*R*(w*da - Delta_omegadelta_w_l);
                xk_l_H_mub = -rhopiR2*R*Nr*h*v*dt;
                xk_l_H_mub_mu0 = -rhopiR2*Nr*h*v*nu0;
                xk_l_H_mub_mu = -rhopiR2*Nr*h*v*Vh;
                xk_l_H_mub_muw = +rhopiR2*Nr*h*v*w;
                xk_l_p = +0.5*rho*(v^2+w^2)*A23*l*(p*l/nu0); % *A23*A12*p/nu0
                
                % My
                xk_m_R_mub = +rhopiR2*R2*v*dr;
                xk_m_T_0 = +rhopiR2*R2*d2e;
                xk_m_T_mu0 = -rhopiR2*R*nu0*de;
                xk_m_T_mu = -rhopiR2*R*Vh*de;
                xk_m_T_muw = +rhopiR2*R*(w*de - Delta_omegadelta_w_m);
                xk_m_H_mub = +rhopiR2*R*Nr*h*u*dt;
                xk_m_H_mub_mu0 = +rhopiR2*Nr*h*u*nu0;
                xk_m_H_mub_mu = +rhopiR2*Nr*h*u*Vh;
                xk_m_H_mub_muw = -rhopiR2*Nr*h*u*w;
                xk_m_q = +0.5*rho*(u^2+w^2)*A13*l*(q*l/nu0); % *A13*A12*q/nu0
                
                % Mz
                xk_n_Q_0 = -rhopiR2*R3*d2r;
                xk_n_Q_mu0 = -rhopiR2*R2*nu0*dr;
                xk_n_Q_mu = -rhopiR2*R2*Vh*dr;
                xk_n_Q_muw = +rhopiR2*R2*(w*dr - Delta_omegadelta_w_n);
                xk_n_Q_muw2 = -rhopiR2*R*Delta_omega_w2_n;
                xk_n_H_mub = -rhopiR2*R*(u*da + v*de); 
                xk_n_H_mub_muw = -rhopiR2*onehalfl2*Nr*(u*p + v*q);
                xk_n_r = +0.5*rho*(u^2+v^2)*((A13+A23)/2)*l*(r*l/nu0); % *A13*A23*r/nu0
            
                % Construct the k'th row of X for each axis
                X{1}(k,:) = [xk_X_H_mub, xk_X_H_mub_mu0, xk_X_H_mub_mu, ...
                             xk_X_H_mub_muw, xk_x_u];
                X{2}(k,:) = [xk_Y_H_mub, xk_Y_H_mub_mu0, xk_Y_H_mub_mu,...
                             xk_Y_H_mub_muw, xk_y_v];
                X{3}(k,:) = [xk_Z_T_0, xk_Z_T_mu0, xk_Z_T_mu, xk_Z_T_mu2,...
                             xk_Z_T_muw, xk_z_w];
                X{4}(k,:) = [xk_l_R_mub, xk_l_T_0, xk_l_T_mu0,xk_l_T_mu,...
                             xk_l_T_muw, xk_l_H_mub, xk_l_H_mub_mu0, xk_l_H_mub_mu,...
                             xk_l_H_mub_muw, xk_l_p];
                X{5}(k,:) = [xk_m_R_mub, xk_m_T_0, xk_m_T_mu0, xk_m_T_mu,...
                             xk_m_T_muw, xk_m_H_mub, xk_m_H_mub_mu0, xk_m_H_mub_mu,...
                             xk_m_H_mub_muw, xk_m_q];
                X{6}(k,:) = [xk_n_Q_0, xk_n_Q_mu0, xk_n_Q_mu, xk_n_Q_muw,...
                             xk_n_Q_muw2, xk_n_H_mub, xk_n_H_mub_muw, xk_n_r];
            
                % Construct the k'th block of H using the k'th rows of X
                %
                % Recall, 
                % theta = [C_H_mub, C_H_mub_mu0, C_H_mub_mu, C_H_mub_muw, C_T_0,...
                %          C_T_mu0, C_T_mu, C_T_mu2, C_T_muw, C_R_mub, C_Q_0, ...
                %          C_Q_mu0, C_Q_mu, C_Q_muw, (C_Q_muw2), C_x_u, C_y_v, ...
                %          C_z_w, C_l_p, C_m_q, C_n_r]'
                %
                H(6*(k-1)+1,[1:4 16]) = X{1}(k,:);
                H(6*(k-1)+2,[1:4 17]) = X{2}(k,:);
                H(6*(k-1)+3,[5:9 18]) = X{3}(k,:);
                H(6*(k-1)+4,[10 5:7 9 1:4 19]) = X{4}(k,:);
                H(6*(k-1)+5,[10 5:7 9 1:4 20]) = X{5}(k,:);
                H(6*(k-1)+6,[11:15 1 4 21]) = X{6}(k,:); 
            
            end
        end % RegressorMatrix

        function detadt = AugmentedDynamics(obj,~,~,~,~,~,~)
            detadt = zeros(0,1);
        end % AugmentedDynamics

        function [dfdv,dfdomega,dfdu,dfdp,dfdeta] = Jacobian(obj,~,vb,omega,Omega,Parameters,Constants)
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
            Vh = sqrt(u^2 + v^2);

            % Angular rates 
            p = omega(1,1);
            q = omega(2,1);
            r = omega(3,1);

            % Virtual actuator states
            delta2 = M*(Omega.^2);
            d2t = delta2(1,1);
            d2a = delta2(2,1);
            d2e = delta2(3,1);
            d2r = delta2(4,1);
            delta = M*Omega;
            dt = delta(1,1);
            da = delta(2,1);
            de = delta(3,1);
            dr = delta(4,1);

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
            C_x_u = C(5);
            C_y_H_mu_b = C(6);
            C_y_H_mu_b_mu = C(7);
            C_y_H_mu_b_mu_w = C(8);
            C_y_H_mu_b_mu_0 = C(9);
            C_y_v = C(10);
            C_z_T_0 = C(11);
            C_z_T_mu = C(12);
            C_z_T_mu_w = C(13);
            C_z_T_mu_0 = C(14);
            C_z_T_mu2 = C(15);
            C_z_w = C(16);
            C_l_T_0 = C(17);
            C_l_T_mu = C(18);
            C_l_T_mu_w = C(19);
            C_l_T_mu_0 = C(20);
            C_l_H_mu_b = C(21);
            C_l_H_mu_b_mu = C(22);
            C_l_H_mu_b_mu_w = C(23);
            C_l_H_mu_b_mu_0 = C(24);
            C_l_R_mu_b = C(25);
            C_l_p = C(26);
            C_m_T_0 = C(27);
            C_m_T_mu = C(28);
            C_m_T_mu_w = C(29);
            C_m_T_mu_0 = C(30);
            C_m_H_mu_b = C(31);
            C_m_H_mu_b_mu = C(32);
            C_m_H_mu_b_mu_w = C(33);
            C_m_H_mu_b_mu_0 = C(34);
            C_m_R_mu_b = C(35);
            C_m_q = C(36);
            C_n_Q_0 = C(37);
            C_n_Q_mu = C(38);
            C_n_Q_mu_w = C(39);
            C_n_Q_mu_0 = C(40);
            C_n_Q_mu_w2 = C(41);
            C_n_H_mu_b = C(42);
            C_n_H_mu_b_mu_w = C(43);
            C_n_r = C(44);

            % The following were computed using symbolic computation on 9/19/2023.
            dfdv = [C_x_u*b*h*rho*u*sign(u)-(C_x_H_mu_b_mu*Nr*rhopiR2*u^2)/Vh-Nr*rhopiR2*(C_x_H_mu_b_mu*Vh+C_x_H_mu_b_mu_0*nu0-C_x_H_mu_b_mu_w*w+C_x_H_mu_b*R*dt)+C_x_u*b*h*rho*u^2*dirac(u),-(C_x_H_mu_b_mu*Nr*rhopiR2*u*v)/Vh,C_x_H_mu_b_mu_w*Nr*rhopiR2*u;
                -(C_y_H_mu_b_mu*Nr*rhopiR2*u*v)/Vh,C_y_v*b*h*rho*v*sign(v)-(C_y_H_mu_b_mu*Nr*rhopiR2*v^2)/Vh-Nr*rhopiR2*(C_y_H_mu_b_mu*Vh+C_y_H_mu_b_mu_0*nu0-C_y_H_mu_b_mu_w*w+C_y_H_mu_b*R*dt)+C_y_v*b*h*rho*v^2*dirac(v),C_y_H_mu_b_mu_w*Nr*rhopiR2*v;
                -rhopiR2*(2*C_z_T_mu2*Nr*u-(C_z_T_mu*Nr*R*dt*u)/Vh),-rhopiR2*(2*C_z_T_mu2*Nr*v-(C_z_T_mu*Nr*R*dt*v)/Vh),C_z_w*b^2*rho*w*sign(w)-C_z_T_mu_w*Nr*R*dt*rhopiR2+C_z_w*b^2*rho*w^2*dirac(w);
                -rhopiR2*((C_l_T_mu*R*da*u)/Vh-C_l_R_mu_b*R^2*dr+(C_l_H_mu_b_mu*Nr*h*u*v)/Vh),(C_l_p*h*p*rho*v*b^3)/nu0-rhopiR2*(C_l_H_mu_b_mu*Nr*h*Vh+C_l_H_mu_b_mu_0*Nr*h*nu0-C_l_H_mu_b_mu_w*Nr*h*w+C_l_H_mu_b*Nr*R*dt*h+(C_l_T_mu*R*da*v)/Vh+(C_l_H_mu_b_mu*Nr*h*v^2)/Vh),(C_l_p*h*p*rho*w*b^3)/nu0+rhopiR2*(C_l_T_mu_w*R*da+C_l_H_mu_b_mu_w*Nr*h*v);
                (C_m_q*h*q*rho*u*b^3)/nu0+rhopiR2*(C_m_H_mu_b_mu*Nr*h*Vh+C_m_H_mu_b_mu_0*Nr*h*nu0-C_m_H_mu_b_mu_w*Nr*h*w+C_m_H_mu_b*Nr*R*dt*h-(C_m_T_mu*R*de*u)/Vh+(C_m_H_mu_b_mu*Nr*h*u^2)/Vh),rhopiR2*(C_m_R_mu_b*R^2*dr-(C_m_T_mu*R*de*v)/Vh+(C_m_H_mu_b_mu*Nr*h*u*v)/Vh),(C_m_q*h*q*rho*w*b^3)/nu0+rhopiR2*(C_m_T_mu_w*R*de-C_m_H_mu_b_mu_w*Nr*h*u);
                (C_n_r*b^3*h*r*rho*u)/nu0-rhopiR2*(C_n_H_mu_b*R*da+C_n_H_mu_b_mu_w*Nr*b^2*p+(C_n_Q_mu*R^2*dr*u)/Vh),(C_n_r*b^3*h*r*rho*v)/nu0-rhopiR2*(C_n_H_mu_b*R*de+C_n_H_mu_b_mu_w*Nr*b^2*q+(C_n_Q_mu*R^2*dr*v)/Vh),C_n_Q_mu_w*R^2*dr*rhopiR2];
            %
            dfdomega = [0,0,0;
                0,0,0;
                C_z_T_mu_w*R*da*rhopiR2,C_z_T_mu_w*R*de*rhopiR2,0;
                (C_l_p*h*rho*(v^2+w^2)*b^3)/(2*nu0)-C_l_T_mu_w*Nr*R*dt*rhopiR2*b^2,C_l_T_mu_w*R*dr*rhopiR2*b^2+J*dr,0;
                C_m_T_mu_w*R*dr*rhopiR2*b^2-J*dr,(C_m_q*h*rho*(u^2+w^2)*b^3)/(2*nu0)-C_m_T_mu_w*Nr*R*dt*rhopiR2*b^2,0;
                -rhopiR2*(-C_n_Q_mu_w*de*R^2+2*C_n_Q_mu_w2*Nr*q*R*b^2+C_n_H_mu_b_mu_w*Nr*u*b^2),-rhopiR2*(-C_n_Q_mu_w*da*R^2+2*C_n_Q_mu_w2*Nr*p*R*b^2+C_n_H_mu_b_mu_w*Nr*v*b^2),(C_n_r*b^3*h*rho*(u^2+v^2))/(2*nu0)];
            %
            dfdu = [0,0,0,0,-C_x_H_mu_b*Nr*R*rhopiR2*u,0,0,0,0;
                0,0,0,0,-C_y_H_mu_b*Nr*R*rhopiR2*v,0,0,0,0;
                -C_z_T_0*Nr*R^2*rhopiR2,0,0,0,rhopiR2*(C_z_T_mu*Nr*R*Vh+C_z_T_mu_0*Nr*R*nu0-C_z_T_mu_w*Nr*R*w),C_z_T_mu_w*R*p*rhopiR2,C_z_T_mu_w*R*q*rhopiR2,0,0;
                0,C_l_T_0*R^2*rhopiR2,0,0,-rhopiR2*(C_l_T_mu_w*Nr*R*p*b^2+C_l_H_mu_b*Nr*R*h*v),-rhopiR2*(C_l_T_mu*R*Vh+C_l_T_mu_0*R*nu0-C_l_T_mu_w*R*w),0,J*q+rhopiR2*(C_l_R_mu_b*u*R^2+C_l_T_mu_w*q*R*b^2),0;
                0,0,C_m_T_0*R^2*rhopiR2,0,rhopiR2*(-C_m_T_mu_w*Nr*R*q*b^2+C_m_H_mu_b*Nr*R*h*u),0,-rhopiR2*(C_m_T_mu*R*Vh+C_m_T_mu_0*R*nu0-C_m_T_mu_w*R*w),rhopiR2*(C_m_R_mu_b*v*R^2+C_m_T_mu_w*p*R*b^2)-J*p,0;
                0,0,0,-C_n_Q_0*R^3*rhopiR2,0,-rhopiR2*(-C_n_Q_mu_w*q*R^2+C_n_H_mu_b*u*R),-rhopiR2*(-C_n_Q_mu_w*p*R^2+C_n_H_mu_b*v*R),-rhopiR2*(C_n_Q_mu*R^2*Vh+C_n_Q_mu_0*R^2*nu0-C_n_Q_mu_w*R^2*w),J];
            %
            dfdp = [-Nr*R*dt*rhopiR2*u,-Nr*rhopiR2*u*Vh,Nr*rhopiR2*u*w,-Nr*nu0*rhopiR2*u,(b*h*rho*u^2*sign(u))/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                0,0,0,0,0,-Nr*R*dt*rhopiR2*v,-Nr*rhopiR2*v*Vh,Nr*rhopiR2*v*w,-Nr*nu0*rhopiR2*v,(b*h*rho*v^2*sign(v))/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                0,0,0,0,0,0,0,0,0,0,-Nr*R^2*d2t*rhopiR2,Nr*R*dt*rhopiR2*Vh,R*rhopiR2*(da*p+de*q-Nr*dt*w),Nr*R*dt*nu0*rhopiR2,-Nr*rhopiR2*(u^2+v^2),(b^2*rho*w^2*sign(w))/2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,R^2*d2a*rhopiR2,-R*da*rhopiR2*Vh,R*rhopiR2*(da*w+b^2*dr*q-Nr*b^2*dt*p),-R*da*nu0*rhopiR2,-Nr*R*dt*h*rhopiR2*v,-Nr*h*rhopiR2*v*Vh,Nr*h*rhopiR2*v*w,-Nr*h*nu0*rhopiR2*v,R^2*dr*rhopiR2*u,(b^3*h*p*rho*(v^2+w^2))/(2*nu0),0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,R^2*d2e*rhopiR2,-R*de*rhopiR2*Vh,R*rhopiR2*(de*w+b^2*dr*p-Nr*b^2*dt*q),-R*de*nu0*rhopiR2,Nr*R*dt*h*rhopiR2*u,Nr*h*rhopiR2*u*Vh,-Nr*h*rhopiR2*u*w,Nr*h*nu0*rhopiR2*u,R^2*dr*rhopiR2*v,(b^3*h*q*rho*(u^2+w^2))/(2*nu0),0,0,0,0,0,0,0,0;
                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-R^3*d2r*rhopiR2,-R^2*dr*rhopiR2*Vh,R^2*rhopiR2*(de*p+da*q+dr*w),-R^2*dr*nu0*rhopiR2,-2*Nr*R*b^2*p*q*rhopiR2,-R*rhopiR2*(da*u+de*v),-Nr*b^2*rhopiR2*(p*u+q*v),(b^3*h*r*rho*(u^2+v^2))/(2*nu0)];
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
                newObj = MultiRotor_BEMT(s.Aircraft,s.Name); 
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