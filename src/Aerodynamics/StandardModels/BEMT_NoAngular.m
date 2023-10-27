classdef BEMT_NoAngular < AerodynamicModel
    %BEMT_NoAngular Summary of this class goes here

    properties
        % Nothing additional.
    end

    methods
    
        function obj = BEMT_NoAngular(Aircraft)
            %AerodynamicModel Construct an instance of this class
            %   Call the AerodynamicModel contstructor.
            obj@AerodynamicModel(Aircraft);
            obj.NumberOfParameters = 17;
        end % BEMT_NoAngular

        function [Force,Moment] = Model(obj,x,Omega,Parameters,Constants)
            %Model

            % Get enviromental constants.
            rho = Constants.rho;
            g = Constants.g;
            
            % Get components of the state vector.
            vb = x(7:9,1);
            omega = x(10:12,1);

            % If C is a parameter estimates object, use the mean. If it is
            % an array, use that.
            if ismatrix(Parameters)
                C = Parameters;
            elseif isobject(Parameters)
                C = Parameters.Mean;
            else
                error('ParamEstimates must be either a ParameterEstimates object or an array.');
            end
            
            % Compute regressor matrix
            RegMat = obj.RegressorMatrix(vb,omega,Omega,rho,g);

            % Compute focres and moments
            % Forces and moments
            XYZLMN = RegMat*C;
            Force = XYZLMN(1:3,1);
            Moment = XYZLMN(4:6,1);
            
        end % Model

        function detadt = AugmentedDynamics(obj,~,~,~,~)
            detadt = [];
        end

        function [dfdx,dfdu,dfdw,dfdp] = Jacobian(obj,x,u,Parameters,Constants)
            %Jacobian

            % Get enviromental constants.
            rho = Constants.rho;
            g = Constants.g;
            
            % Get components of the state vector.
            vb = x(7:9,1);
            omega = x(10:12,1);
            
            % For a linearly parameterized model, dfdp is the regressor
            % matrix padded with zeros.
            dFMdp = obj.RegressorMatrix(vb,omega,u,rho,g);
            dfdp = [zeros(6,size(dFMdp,2));dFMdp];

            % Compute the rigid body jacobian.
            J = obj.Aircraft.RigidBodyJacobian(x,g);

            % Compute the derivatives of the aerodynamic model with respect
            % to states and inputs.
            %
            % For now, comput this numerically. TODO: The derivative is not
            % bad to analytically compute. I have already done this in
            % Mathematica...
            %
            epsilon = 1e-8;
            epsilon_inv = 1/epsilon;
            dFMdx = zeros(6,12);
            for ii = 1:12
                xplus = x;
                xminus = x;
                xplus(ii,1) =  x(ii,1) + epsilon;
                xminus(ii,1) =  x(ii,1) - epsilon;
                [Fplus,Mplus] = obj.Model(xplus,u,Parameters,Constants);
                [Fminus,Mminus] = obj.Model(xminus,u,Parameters,Constants);
                dFMdx(:,ii) = 0.5*epsilon_inv*[Fplus-Fminus;Mplus-Mminus];
            end
            %
            nu = obj.Aircraft.NumberOfRotors;
            dFMdu = zeros(6,obj.Aircraft.NumberOfRotors);
            for ii = 1:nu
                uplus = u;
                uminus = u;
                uplus(ii,1) =  u(ii,1) + epsilon;
                uminus(ii,1) =  u(ii,1) - epsilon;
                [Fplus,Mplus] = obj.Model(x,uplus,Parameters,Constants);
                [Fminus,Mminus] = obj.Model(x,uminus,Parameters,Constants);
                dFMdu(:,ii) = 0.5*epsilon_inv*[Fplus-Fminus;Mplus-Mminus];
            end

            % Construct the Jacobains of the aircraft dynamics
            dfdx = J + [zeros(6,12);dFMdx];
            dfdu = [zeros(6,4);dFMdu];
            dfdw = [zeros(6,3);dFMdx(:,7:9)];
            
        end % Model
      
        function [X,TexNames] = RegressorMatrix(obj,vb,omega,Omega,rho,g)
            %RegressorMatrix

            % aircraft mass
            m = obj.Aircraft.Mass;
            
            % aircraft geometry
            R = obj.Aircraft.RotorDiameter/2;
            M_mix = obj.Aircraft.MixingMatrix;
            N_rot = obj.Aircraft.NumberOfRotors;
            b = obj.Aircraft.ArmLength*sqrt(2)/2;
            h = obj.Aircraft.HubHeight;
            
            % emperical terms
            nu0 = sqrt(m*g./(2*rho*pi*R^2));
            
            % explanatory variables
            u = vb(1,1);
            v = vb(2,1);
            w = vb(3,1);
            V = sqrt(u.^2+v.^2+w.^2);
            Vh = sqrt(u.^2+v.^2);
            % p = omega(1,1);
            % q = omega(2,1);
            % r = omega(3,1);
            delta = M_mix*Omega.^2;
            deltaSqrt = M_mix*Omega;
            dt = delta(1,1);
            da = delta(2,1);
            de = delta(3,1);
            dr = delta(4,1);
            sqdt = deltaSqrt(1,1);
            sqda = deltaSqrt(2,1);
            sqde = deltaSqrt(3,1);
            sqdr = deltaSqrt(4,1);
            
            % shared terms
            rho_pi_R2 = rho.*pi.*R.^2;
            qbar_b2 = 0.5.*rho.*V.^2.*b.^2;
            
            % candidate regressors
            X_H_mub = -rho_pi_R2.*R.*u.*sqdt;
            X_H_mubmu = -rho_pi_R2.*N_rot.*u.*Vh;
            X_H_mubmuw = -rho_pi_R2.*N_rot.*u.*w;
            X_H_mubmu0 = -rho_pi_R2.*N_rot.*u.*nu0;
            X_u = qbar_b2.*u;
            Y_H_mub = -rho_pi_R2.*R.*v.*sqdt;
            Y_H_mubmu = -rho_pi_R2.*N_rot.*v.*Vh;
            Y_H_mubmuw = -rho_pi_R2.*N_rot.*v.*w;
            Y_H_mubmu0 = -rho_pi_R2.*N_rot.*v.*nu0;
            Y_v = qbar_b2.*v;
            Z_T_0 = -rho_pi_R2.*R^2.*dt;
            Z_T_mu0 = +rho_pi_R2.*R.*nu0.*sqdt;
            Z_T_mu = +rho_pi_R2.*R.*Vh.*sqdt;
            Z_T_muw = +rho_pi_R2.*R.*w.*sqdt;
            Z_T_mu2 = -rho_pi_R2.*N_rot.*Vh.^2;
            Z_w = qbar_b2.*w;
            L_T_0 = +rho_pi_R2.*b.*R.^2.*da;
            L_T_mu0 = -rho_pi_R2.*b.*R.*nu0.*sqda;
            L_T_mu = -rho_pi_R2.*b.*R.*Vh.*sqda;
            L_T_muw = -rho_pi_R2.*b.*R.*w.*sqda;
            L_H_mub = -rho_pi_R2.*h.*R.*v.*sqdt;
            L_H_mubmu0 = -rho_pi_R2.*N_rot.*h.*v.*nu0;
            L_H_mubmu = -rho_pi_R2.*N_rot.*h.*v.*Vh;
            L_H_mubmuw = -rho_pi_R2.*N_rot.*h.*v.*w;
            L_R_mub = -rho_pi_R2.*R.^2.*u.*sqdr;
            M_T_0 = +rho_pi_R2.*b.*R.^2.*de;
            M_T_mu0 = -rho_pi_R2.*b.*R.*nu0.*sqde;
            M_T_mu = -rho_pi_R2.*b.*R.*Vh.*sqde;
            M_T_muw = -rho_pi_R2.*b.*R.*w.*sqde;
            M_H_mub = +rho_pi_R2.*h.*R.*u.*sqdt;
            M_H_mubmu0 = +rho_pi_R2.*N_rot.*h.*u.*nu0;
            M_H_mubmu = +rho_pi_R2.*N_rot.*h.*u.*Vh;
            M_H_mubmuw = +rho_pi_R2.*N_rot.*h.*u.*w;
            M_R_mub = -rho_pi_R2.*R.^2.*v.*sqdr;
            N_Q_0 = -rho_pi_R2.*R.^3.*dr;
            N_Q_mu0 = -rho_pi_R2.*R.^2.*nu0.*sqdr;
            N_Q_mu = -rho_pi_R2.*R.^2.*Vh.*sqdr;
            N_Q_muw = -rho_pi_R2.*R.^2.*w.*sqdr;
            N_H_mub = -rho_pi_R2.*b.*R.*(u.*sqda+v.*sqde);
            
            % Regressor matrix for linearly parameterized model
            X = [zeros(1,5) X_H_mub X_H_mubmu0 X_H_mubmu X_H_mubmuw zeros(1,5) X_u zeros(1,2);
                 zeros(1,5) Y_H_mub Y_H_mubmu0 Y_H_mubmu Y_H_mubmuw zeros(1,6) Y_v zeros(1,1);
                 Z_T_0 Z_T_mu0 Z_T_mu Z_T_muw Z_T_mu2 zeros(1,11) Z_w;
                 L_T_0 L_T_mu0 L_T_mu L_T_muw zeros(1,1) L_H_mub L_H_mubmu0 L_H_mubmu L_H_mubmuw zeros(1,4) L_R_mub zeros(1,3);
                 M_T_0 M_T_mu0 M_T_mu M_T_muw zeros(1,1) M_H_mub M_H_mubmu0 M_H_mubmu M_H_mubmuw zeros(1,4) M_R_mub zeros(1,3);
                 zeros(1,5) N_H_mub zeros(1,3) N_Q_0 N_Q_mu0 N_Q_mu N_Q_muw zeros(1,4)];

            % parameter names
            TexNames = {'$C_{T_0}$';
                        '$C_{T_{\mu_0}}$';
                        '$C_{T_\mu}$';
                        '$C_{T_{\mu_w}}$';
                        '$C_{T_{\mu^2}}$';
                        '$C_{H_{\mu_b}}$';
                        '$C_{H_{\mu_b,\mu_0}}$';
                        '$C_{H_{\mu_b,\mu}}$';
                        '$C_{H_{\mu_b,\mu_w}}$';
                        '$C_{Q_0}$';
                        '$C_{Q_{\mu_0}}$';
                        '$C_{Q_\mu}$';
                        '$C_{Q_{\mu_w}}$';
                        '$C_{R_{\mu_b}}$';
                        '$C_{x_u}$';
                        '$C_{y_v}$';
                        '$C_{z_w}$'};

        end % RegressorMatrix

    end % public methods

    methods(Static)

        function obj = loadobj(s)
            %loadobj
            %
            % If any of the saved properties no longer are defined in the
            % class definition, create a dynamic property and give a
            % warning.
            if isstruct(s)
                newObj = BEMT_NoAngular(s.Aircraft); 
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