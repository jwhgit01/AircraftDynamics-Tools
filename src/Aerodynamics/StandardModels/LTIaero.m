classdef LTIaero < AerodynamicModel
%LTIaero

properties
    % Nothing additional
end

methods

    function obj = LTIaero(Aircraft,Name)
        %LTIaero Construct an instance of this class
        %   Call the AerodynamicModel contstructor.
        obj@AerodynamicModel(Aircraft,Name);

        % Set the number of dynamic states and parameters
        obj.NumberOfStates = 0;
        obj.NumberOfParameters = 66;

        % check to see if the CSV already exists
        AircraftName = obj.Aircraft.Name;
        AircraftPath = obj.Aircraft.WorkspaceDir;
        suffix = 'Parameters_LTIaero.csv';
        if 0 ~= exist([AircraftPath '/AircraftData/' AircraftName '_' suffix],'file')
            warning(['A LTIaero model already exists for ' AircraftName]);
            return
        end
        
        % Copy the template CSV
        temp = which(['AircraftDynamcis-Tools/lib/' suffix]);
        [libpath,~,~] = fileparts(temp);
        copyfile([libpath filesep suffix],[AircraftPath '/AircraftData/' AircraftName '_' suffix]);

    end % LTIaero

    function [F,M] = Model(obj,~,vb_r,omega_r,delta,Parameters,~)
        %Model

        % Get parameters and construct LTI aerodynamic matrices
        if isobject(Parameters)
            C = Parameters.Mean;
        else
            C = Parameters;
        end
        C_0 = C(1:3,1);
        D_0 = C(4:6,1);
        C_v = reshape(C(7:15,1),3,3).';
        D_v = reshape(C(16:24,1),3,3).';
        C_omega = reshape(C(25:33,1),3,3).';
        D_omega = reshape(C(34:42,1),3,3).';
        C_delta = reshape(C(43:54,1),4,3).';
        D_delta = reshape(C(55:66,1),4,3).';

        % Compute forces and moments
        F = -C_0 + C_v*vb_r + C_omega*omega_r + C_delta*delta;
        M = -D_0 + D_v*vb_r + D_omega*omega_r + D_delta*delta;
        
    end % Model

    function detadt = AugmentedDynamics(obj,~,~,~,~,~,~)
        detadt = zeros(0,1);
    end

    function [dfdv,dfdomega,dfddelta,dfdp,dfdeta] = Jacobian(obj,~,vb_r,omega_r,delta,Parameters,~)
        %Jacobian

        % Get parameters and construct LTI aerodynamic matrices
        if isobject(Parameters)
            C = Parameters.Mean;
        else
            C = Parameters;
        end
        C_v = reshape(C(7:15,1),3,3).';
        D_v = reshape(C(16:24,1),3,3).';
        C_omega = reshape(C(25:33,1),3,3).';
        D_omega = reshape(C(34:42,1),3,3).';
        C_delta = reshape(C(43:54,1),4,3).';
        D_delta = reshape(C(55:66,1),4,3).';
            
        % Construct Jacobians
        dfdv = [C_v; D_v];
        dfdomega = [C_omega; D_omega];
        dfddelta = [C_delta; D_delta];
        dfdp_0 = eye(6);
        dfdp_v = blkdiag(vb_r.',vb_r.',vb_r.',vb_r.',vb_r.',vb_r.');
        dfdp_omega = blkdiag(omega_r.',omega_r.',omega_r.',omega_r.',omega_r.',omega_r.');
        dfdp_delta = blkdiag(delta.',delta.',delta.',delta.',delta.',delta.');
        dfdp = [dfdp_0, dfdp_v, dfdp_omega, dfdp_delta];
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
            newObj = LTIaero(s.Aircraft,s.Name); 
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