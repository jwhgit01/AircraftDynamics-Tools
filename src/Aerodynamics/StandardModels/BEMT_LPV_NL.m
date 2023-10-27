classdef BEMT_LPV_NL < AerodynamicModel
    %BEMT_LPV_NL
    %
    % This model is the nonlinear model which upon linearization, becomes
    % BEMT_LPV.m

    properties
        % Nothing additional.
    end

    methods
    
        function obj = BEMT_LPV_NL(Aircraft)
            %AerodynamicModel Construct an instance of this class
            %   Call the AerodynamicModel contstructor.
            obj@AerodynamicModel(Aircraft);
        end % BEMT_NoAngular

        function [Force,Moment] = Model(obj,x,Omega,Parameters,Constants)
            %Model

            % Get enviromental constants.
            rho = Constants.rho;
            
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
            [HX,HY,HZ,HL,HM,HN,~] = obj.RegressorMatrix(vb,omega,Omega,rho);

            % Compute forces and moments
            n1 = size(HX,2);
            n2 = size(HY,2) + n1;
            n3 = size(HZ,2) + n2;
            n4 = size(HL,2) + n3;
            n5 = size(HM,2) + n4;
            np = size(HN,2) + n5;
            obj.NumberOfParameters = np;
            X = HX*C(1:n1);
            Y = HY*C(n1+1:n2);
            Z = HZ*C(n2+1:n3);
            L = HL*C(n3+1:n4);
            M = HM*C(n4+1:n5);
            N = HN*C(n5+1:np);
            Force = [X;Y;Z];
            Moment = [L;M;N];
            
        end % Model

        function detadt = AugmentedDynamics(obj,~,~,~,~)
            detadt = [];
        end

        function [dfdx,dfdu,dfdw,dfdp] = Jacobian(obj,x,u,Parameters,Constants)
            %Jacobian

            % Get enviromental constants.
            rho = Constants.rho;
            g = Constants.g;
            
            % Get components of the state vector and virtual actuators.
            vb = x(7:9,1);
            omega = x(10:12,1);
            M_mix = obj.Aircraft.MixingMatrix;
            deltaSqrt = M_mix*u;
            sqdt = deltaSqrt(1,1);
            sqda = deltaSqrt(2,1);
            sqde = deltaSqrt(3,1);
            sqdr = deltaSqrt(4,1);
            
            % For a linearly parameterized model, dfdp is the regressor
            % matrix padded with zeros.
            [HX,HY,HZ,HL,HM,HN,~] = obj.RegressorMatrix(vb,omega,u,rho);
            dFMdp = blkdiag(HX,HY,HZ,HL,HM,HN);
            np = size(dFMdp,2);
            dfdp = [zeros(6,np);dFMdp];

            % Compute the rigid body jacobian.
            J = obj.Aircraft.RigidBodyJacobian(x,g);

            % Parameters
            mc = metaclass(Parameters);
            if strcmp(mc.Name,'ParameterEstimates')
                C = Parameters.Mean;
            else
                C = Parameters;
            end

            % Compute the derivatives of the aerodynamic model with respect
            % to states and inputs.
            dFMdx = zeros(6,12,'like',C);
            dFMdx(1,7) = C(1)*sqdt + C(2)*vb(3,1);
            dFMdx(1,9) = C(2)*vb(1,1);
            dFMdx(2,8) = C(3)*sqdt + C(4)*vb(3,1);
            dFMdx(2,9) = C(4)*vb(2,1);
            dFMdx(3,7) = C(6)*vb(1,1)*sqdt;
            dFMdx(3,8) = C(7)*vb(2,1)*sqdt;
            dFMdx(3,9) = C(8)*sqdt;
            dFMdx(4,7) = C(13)*sqdr;
            dFMdx(4,8) = C(11)*sqdt + C(12)*vb(3,1);
            dFMdx(4,9) = C(10)*sqda + C(12)*vb(2,1);
            dFMdx(4,10) = C(14);
            dFMdx(5,7) = C(17)*sqdt + C(18)*vb(3,1);
            dFMdx(5,8) = C(19)*sqdr;
            dFMdx(5,9) = C(16)*sqde + C(18)*vb(1,1);
            dFMdx(5,11) = C(20);
            dFMdx(6,7) = C(22)*vb(1,1)*sqdr + C(25)*sqda;
            dFMdx(6,8) = C(23)*vb(2,1)*sqdr + C(26)*sqde;
            dFMdx(6,9) = C(24)*sqdr;
            dFMdx(6,12) = C(27);
            %
            dFMdu = zeros(6,4,'like',C);
            dFMdu(1,1) = C(1)*vb(1,1);
            dFMdu(2,1) = C(3)*vb(2,1);
            dFMdu(3,1) = C(5) + C(8)*vb(3,1);
            dFMdu(4,1) = C(11)*vb(2,1);
            dFMdu(4,2) = C(9) + C(10)*vb(3,1);
            dFMdu(4,4) = C(13)*vb(1,1);
            dFMdu(5,1) = C(17)*vb(1,1);
            dFMdu(5,3) = C(15) + C(16)*vb(3,1);
            dFMdu(5,4) = C(19)*vb(2,1);
            dFMdu(6,2) = C(25)*vb(1,1);
            dFMdu(6,3) = C(26)*vb(2,1);
            dFMdu(6,4) = C(21) + C(24)*vb(3,1);

            % Construct the Jacobains of the aircraft dynamics
            dfdx = J + [zeros(6,12);dFMdx];
            dfdu = [zeros(6,4);dFMdu];
            dfdw = [zeros(6,6);dFMdx(:,7:12)];
            
        end % Jacobian

        function [y,H,pidx] = StateSpaceEquationError(obj,dx,du,x0,u0,accel,Constants)
            %StateSpaceEquationError   

            % Which elements of A and B contain parameters from the 
            % hypothesized LPV model structure?
            pIdxA0 = [7 7 8 8 9 9 9 10 10 10 10 11 11 11 11 12 12 12 12;
                      7 9 8 9 7 8 9  7  8  9 10  7  8  9 11  7  8  9 12];
            pIdxB0 = [7 8 9 10 10 10 11 11 11 12 12 12;
                      1 1 1  1  2  4  1  3  4  2  3  4];

            % Check to see if any of these parameters vanish at the current
            % linearization point.
            np = obj.NumberOfParameters;
            [dfdx,dfdu,~,~] = Jacobian(obj,x0,u0,ones(np,1),Constants);
            pIdxA = pIdxA0; 
            jj = 0;
            for ii = 1:size(pIdxA0,2)
                if dfdx(pIdxA0(:,ii)) ~= 0
                    jj = jj + 1;
                    pIdxA(:,jj) = pIdxA0(:,ii);
                end
            end
            pIdxB = pIdxB0;
            jj = 0;
            for ii = 1:size(pIdxB0,2)
                if dfdu(pIdxB0(:,ii)) ~= 0
                    jj = jj + 1;
                    pIdxB(:,jj) = pIdxB0(:,ii);
                end
            end
            
            % Loop through the parameter indices and construct regressors.
            % Also, create array of indices that indicate the element of
            % the A/B matrix the parameter represents.
            H = cell(6,1);
            y = zeros(size(dx,1),6);
            pidx = [];
            for ii = 1:6
                % Logical indices of the state/input vector that are
                % regressors for the ith axis.
                xidx = pIdxA(1,:)==6+ii;
                uidx = pIdxB(1,:)==6+ii;
                % Indices of the A and B matrices that correspond to these
                % regressors.
                pidxAi = pIdxA(:,xidx);
                pidxBi = pIdxB(:,uidx);
                % The number of parameters in the A and B matrices for this
                % axis.
                npAi = size(pidxAi,2);
                npBi = size(pidxBi,2);
                % An array of indices indicating each parameters
                % corresponding matrix and indices.
                pidxi = [ones(1,npAi), 2*ones(1,npBi)
                               pidxAi,        pidxBi];
                pidx = [pidx pidxi];
                % Regressor matrix for the ith axis.
                HidxA = pIdxA(2,xidx);
                HidxB = pIdxB(2,uidx);
                H{ii} = [dx(:,HidxA), du(:,HidxB)];
                % Measurement
                yidxA = setxor(HidxA,1:12);
                yidxB = setxor(HidxB,1:4);
                y(:,ii) = accel(:,ii)-(dfdx(ii,yidxA)*(dx(:,yidxA).')+dfdu(ii,yidxB)*(du(:,yidxB).')).';
            end
            
        end % StateSpaceRegressorMatrix

        function [A,B,C,D,K] = GreyBox(obj,Parameters,~,Constants)
            %GreyBox   
            %
            % Has the form:
            %   [A,B,C,D,K] = myfunc(par1,par2,...,parN,Ts,aux1,aux2,...)

            % Compute Jacobians
            [dfdx,dfdu,dfdw,~] = Jacobian(obj,Constants.x0,Constants.u0,Parameters,Constants);

            % Construct linear system matrices
            A = dfdx;
            B = dfdu;
            C = eye(12);
            D = zeros(12,4);
            K = dfdw;
            
        end % GreyBox
      
        function [HX,HY,HZ,HL,HM,HN,TexNames] = RegressorMatrix(obj,vb,omega,Omega,rho)
            %RegressorMatrix
            
            % aircraft geometry
            R = obj.Aircraft.RotorDiameter/2;
            M_mix = obj.Aircraft.MixingMatrix;
            
            % explanatory variables
            u = vb(1,1);
            v = vb(2,1);
            w = vb(3,1);
            p = omega(1,1);
            q = omega(2,1);
            r = omega(3,1);
            deltaSqrt = M_mix*Omega;
            sqdt = deltaSqrt(1,1);
            sqda = deltaSqrt(2,1);
            sqde = deltaSqrt(3,1);
            sqdr = deltaSqrt(4,1);
            
            % shared terms
            rho_pi_R2 = rho.*pi.*R.^2;
            
            % candidate regressors
            X_mub = rho_pi_R2.*u.*sqdt;
            X_mubmuw = rho_pi_R2.*u.*w;
            Y_mub = rho_pi_R2.*v.*sqdt;
            Y_mubmuw = rho_pi_R2.*v.*w;
            Z_T = rho_pi_R2.*sqdt;
            Z_muu = 0.5*rho_pi_R2.*u.^2.*sqdt;
            Z_muv = 0.5*rho_pi_R2.*v.^2.*sqdt;
            Z_muw = rho_pi_R2.*w.*sqdt;
            L_T = rho_pi_R2.*sqda;
            L_muw = rho_pi_R2.*w.*sqda;
            L_mub = rho_pi_R2.*v.*sqdt;
            L_mubmuw = rho_pi_R2.*v.*w;
            L_R = rho_pi_R2.*u.*sqdr;
            L_p = rho_pi_R2.*p;
            M_T = rho_pi_R2.*sqde;
            M_muw = rho_pi_R2.*w.*sqde;
            M_mub = rho_pi_R2.*u.*sqdt;
            M_mubmuw = rho_pi_R2.*u.*w;
            M_R = rho_pi_R2.*v.*sqdr;
            M_q = rho_pi_R2.*q;
            N_Q = rho_pi_R2.*sqdr;
            N_muu = 0.5*rho_pi_R2.*u.^2.*sqdr;
            N_muv = 0.5*rho_pi_R2.*v.^2.*sqdr;
            N_muw = rho_pi_R2.*w.*sqdr;
            N_Hu = rho_pi_R2.*u.*sqda;
            N_Hv = rho_pi_R2.*v.*sqde;
            N_r = rho_pi_R2.*r;
            
            % Regressor matrices for linearly parameterized model
            HX = [X_mub, X_mubmuw];
            HY = [Y_mub, Y_mubmuw];
            HZ = [Z_T, Z_muu, Z_muv Z_muw];
            HL = [L_T, L_muw, L_mub, L_mubmuw, L_R, L_p];
            HM = [M_T, M_muw, M_mub, M_mubmuw, M_R, M_q];
            HN = [N_Q, N_muu, N_muv, N_muw, N_Hu, N_Hv, N_r];

            % parameter names
            TexNames.X = {'$C_{X_{\mu_b}}$';'$C_{X_{\mu_b,\mu_w}}$'};
            TexNames.Y = {'$C_{Y_{\mu_b}}$';'$C_{Y_{\mu_b,\mu_w}}$'};
            TexNames.Z = {'$C_{Z_{T}}$';'$C_{Z_{\mu_u}}$';'$C_{Z_{\mu_v}}$';'$C_{Z_{\mu_w}}$'};
            TexNames.L = {'$C_{l_{T}}$';'$C_{l_{\mu_w}}$';'$C_{l_{\mu_b}}$';'$C_{l_{\mu_b,\mu_w}}$';'$C_{l_{R}}$';'$C_{l_{p}}$'};
            TexNames.M = {'$C_{m_{T}}$';'$C_{m_{\mu_w}}$';'$C_{m_{\mu_b}}$';'$C_{m_{\mu_b,\mu_w}}$';'$C_{m_{R}}$';'$C_{m_{q}}$'};
            TexNames.N = {'$C_{n_{Q}}$';'$C_{n_{\mu_u}}$';'$C_{n_{\mu_v}}$';'$C_{n_{\mu_w}}$';'$C_{n_{H_u}}$';'$C_{n_{H_v}}$';'$C_{n_{r}}$'};

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
                newObj = BEMT_LPV_NL(s.Aircraft); 
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