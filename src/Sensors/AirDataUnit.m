classdef AirDataUnit < Sensor
    %AirDataUnit Summary of this class goes here

    properties
        Position % geometric center of vanes/pitot tube (body frame)
        Bias % struct with fields {PitotTube, BetaVane, AlphaVane}
        DynamicModel
    end

    methods
        function obj = AirDataUnit(Name,Position)
            %AirDataUnit Construct an instance of this class
            obj.Name = Name;
            obj.Position = Position;
        end % AirDataUnit

        function y = Model(obj,x,~,~,~)
            %Model y = [V_ind;BetaVane;AlphaVane]

            vb = x(7:9,1);
            omega = x(10:12,1);
            vs = vb + cross(omega,obj.Position);
            alpha_s = atan2(vs(3),vs(2));
            mu_s = atan2(vs(2),vs(1));
            V_s = norm(vs);
            V_ind = V_s;
            BetaVane = mu_s;
            AlphaVane = alpha_s;
            y = [V_ind; BetaVane; AlphaVane];

            if ~isempty(obj.Bias)
                y = y + obj.Bias;
            end

        end % Model

        function [dhdx,dhdu,dhdp] = Jacobian(obj,x,u,Parameters,Constants)
            %Jacobian Measurement model Jacobians with respect to the state
            % vector, x = [NED; Theta; vb; omega; eta], and the parameter
            % vector, p. Here, eta are any additional states used in the
            % aerodynamic model.

            dhdx = [];
            dhdu = [];
            dhdp = [];

        end % Jacobian

        function vbr = InverseModel(obj,V_ind,BetaVane,AlphaVane,omega)
            %InverseModel
            vbr = zeros(length(V_ind),3);
            for ii = 1:length(V_ind)
                V_s = V_ind(ii);
                mu_s = BetaVane(ii);
                alpha_s = AlphaVane(ii);
                beta_s = atan(tan(mu_s).*cos(alpha_s));
                vs = [V_s.*cos(alpha_s).*cos(beta_s);
                      V_s.*sin(beta_s);
                      V_s.*sin(alpha_s).*cos(beta_s)];
                vbr(ii,:) = (vs - cross(omega(ii,:).',obj.Position)).';
            end
        end % InverseModel
      
    end % public methods

    methods(Static)
        
        function obj = loadobj(s)
            %loadobj
            %
            % If any of the saved properties no longer are defined in the
            % class definition, create a dynamic property and give a
            % warning.
            if isstruct(s)
                newObj = AirDataUnit(s.Name,s.Position); 
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