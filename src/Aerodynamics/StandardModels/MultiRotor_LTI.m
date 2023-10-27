classdef MultiRotor_LTI < AerodynamicModel
    %MultiRotor_LTI

    properties
        
        % The x0 and u0 properties define the point about which this linear
        % model is deifined, where x0 is the nominal state vector and u0 is
        % the nominal virtual actuator input vector.
        x0
        u0

        % The ParameterIndices property specifies the element of the A or B
        % matrix each parameter represents. The first row of
        % ParameterIndices specifies whether the parameter defines an
        % element of the A (=1) or B matrix (=2). The seconds and third
        % rows define the row and column indices, respectively.
        ParameterIndices

    end

    methods
    
        function obj = MultiRotor_LTI(Aircraft,Name)
            %MultiRotor_LTI Construct an instance of this class
            %   Call the AerodynamicModel contstructor.
            obj@AerodynamicModel(Aircraft,Name);
        end % MultiRotor_LTI

        function [Force,Moment] = Model(obj,x,u,Parameters,Constants)
            %Model
            %
            % This model is only valid for the case when (x0,u0) is such
            % that dv/dt = domega/dt = 0.

            % Get constants.
            m = obj.Aircraft.Mass;
            I = obj.Aircraft.MomentOfInertia;
            g = Constants.g;

            % Get states.
            vb = x(7:9,1);
            omega = x(10:12,1);
            
            % Compute Jacobians.
            [A,B,~,~] = Jacobian(obj,x,u,Parameters,Constants);

            % State perturbations.
            dx = x - obj.x0;
            du = u - obj.u0;
            
            % Subtract the contribution from the rigid body dynamics from
            % the linear state equation to get expressions for the forces
            % and moments.
            e3 = [0;0;1];
            phi = x(4,1);
            theta = x(5,1);
            psi = x(6,1);
            R1 = [1,0,0;0,cos(phi),-sin(phi);0,sin(phi),cos(phi)];
            R2 = [cos(theta),0,sin(theta);0,1,0;-sin(theta),0,cos(theta)];
            R3 = [cos(psi),-sin(psi),0;sin(psi),cos(psi),0;0,0,1];
            R_BI = (R3(psi)*R2(theta)*R1(phi)).';
            Force = cross(omega,m*vb) - m*g*R_BI*e3 + m*A(7:9,:)*dx + m*B(7:9,:)*du;
            Moment = cross(omega,I*omega) + I*(A(10:12,:)*dx + B(10:12,:)*du);
            
        end % Model

        function detadt = AugmentedDynamics(obj,~,~,~,~)
            detadt = [];
        end

        function [dfdx,dfdu,dfdw,dfdp] = Jacobian(obj,x,u,Parameters,Constants)
            %Jacobian

            % Get enviromental and aircraft constants.
            rho = Constants.rho;
            g = Constants.g;
            R = obj.Aircraft.RotorDiameter/2;
            rhopiR2 = rho*pi*R^2;

            % Parameters
            mc = metaclass(Parameters);
            if strcmp(mc.Name,'ParameterEstimates')
                C = Parameters.Mean;
            else
                C = Parameters;
            end
            
            % Populate the A and B matrices based on the ParameterIndices
            % property. Remember, parameters are scaled by rho*pi*R^2.
            pidx = obj.ParameterIndices;
            np = size(pidx,2);
            dfdx = obj.Aircraft.RigidBodyJacobian(obj.x0,g);
            dfdu = zeros(12,4);
            if strcmp(mc.Name,'umat')
                dfdx = umat(dfdx);
                dfdu = umat(dfdu);
            end
            for pp = 1:np
                ii = pidx(2,pp);
                jj = pidx(3,pp);
                if pidx(1,pp) == 1
                    dfdx(ii,jj) = rhopiR2*C(pp);
                elseif pidx(1,pp) == 2
                    dfdu(ii,jj) = rhopiR2*C(pp);
                else
                    error('Invalid first row of ParameterIndices property')
                end
            end
            
            % The system disturbances are assumed to be forces and moments.
            dfdw = [zeros(6,6);eye(6)];

            % The jacobian of the dynamics with repect to the parameters is
            % constructed using the ParameterIndices property as well.
            dfdp = zeros(12,np);
            for pp = 1:np
                ii = pidx(2,pp);
                jj = pidx(3,pp);
                if pidx(1,pp) == 1
                    dfdp(ii,pp) = rhopiR2*x(jj);
                elseif pidx(1,pp) == 2
                    dfdp(ii,pp) = rhopiR2*u(jj);
                else
                    error('Invalid first row of ParameterIndices property')
                end
            end
            
        end % Jacobian

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

        function [y,H] = EquationError(obj,x,u,accel,Constants)
            %EquationError  

            % Get enviromental and aircraft constants.
            rho = Constants.rho;
            R = obj.Aircraft.RotorDiameter/2;
            rhopiR2 = rho*pi*R^2;

            % State and input perturbations
            dx = x - obj.x0.';
            du = u - obj.u0.';

            % Parameter indices and number of parameters.
            pidx = obj.ParameterIndices;
            np = size(pidx,2);

            % Compute Jacobians (ignore elements with parameters)
            [A0,B0,~,~] = Jacobian(obj,x,u,zeros(np,1),Constants);

            % Construct parameter index sets.
            IA = pidx(2:3,pidx(1,:)==1);
            IB = pidx(2:3,pidx(1,:)==2);
            [row,col] = find(A0);
            IAbar = [row col]';
            [row,col] = find(B0);
            IBbar = [row col]';
            
            % Loop through state equations and construct regressors and
            % measurements.
            H = cell(6,1);
            y = zeros(size(x,1),6);
            for ii = 1:6

                % Use compliment index sets to construct measurements.
                jjA = IAbar(2,IAbar(1,:)==6+ii);
                jjB = IBbar(2,IBbar(1,:)==6+ii);
                y(:,ii) = accel(:,ii) - (A0(6+ii,jjA)*(dx(:,jjA).')).' - (B0(6+ii,jjB)*(du(:,jjB).')).';

                % Use index sets to construct regressor matrices.
                jjA = IA(2,IA(1,:)==6+ii);
                jjB = IB(2,IB(1,:)==6+ii);
                H{ii} = rhopiR2*[dx(:,jjA), du(:,jjB)];

            end
            
        end % EquationError

    end % public methods

    methods(Static)
        
        function obj = loadobj(s)
            %loadobj
            %
            % If any of the saved properties no longer are defined in the
            % class definition, create a dynamic property and give a
            % warning.
            if isstruct(s)
                newObj = BEMT_LTI(s.Aircraft); 
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