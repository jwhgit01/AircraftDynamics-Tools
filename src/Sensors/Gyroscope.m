classdef Gyroscope < Sensor
    %Gyroscope Summary of this class goes here
    %
    % Assumptions:
    %   - The gyro bias is constant.

    properties
        Bias
    end

    methods
        function obj = Gyroscope(Name)
            %Gyroscope Construct an instance of this class
            obj.Name = Name;
        end % Gyroscope

        function y = Model(obj,x,~,~,~)
            %Model Measurement model for a gyroscope. The output, y,
            % is y = [omegax; omegay; omegaz] + bias.
            
            % Angular velocity
            y = x(10:12,1); 

            % If there is bias, add it.
            if ~isempty(obj.Bias)
                y = y + obj.Bias;
            end

        end

        function [dhdx,dhdu,dhdp] = Jacobian(obj,x,u,Parameters,~)
            %Jacobian Measurement model Jacobian with respect to the state
            %vector, x = [NED; Theta; vb; omega; eta].

            % Get necessary dimensions.
            nx = size(x,1);
            nu = size(u,1);
            if isobject(Parameters)
                np = size(Parameters.Mean,1);
            else
                np = size(Parameters,1);
            end

            % Jacobians are constant.
            dhdx = [zeros(3,9), eye(3), zeros(3,nx-12)];
            dhdu = zeros(3,nu);
            dhdp = zeros(3,np);
            
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
                newObj = Gyroscope(s.Name); 
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