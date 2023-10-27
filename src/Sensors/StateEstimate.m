classdef StateEstimate < Sensor
    %StateEstimate Summary of this class goes here
    %

    properties
        Index
    end

    methods
        function obj = StateEstimate(Name,Index)
            %StateEstimate Construct an instance of this class
            obj.Name = Name;
            obj.Index = Index;
        end % StateEstimate

        function y = Model(obj,x,~,~,~)
            %Model Measurement model.
            y = x(obj.Index,1);
        end

        function [dhdx,dhdu,dhdp] = Jacobian(obj,x,u,Parameters,~)
            %Jacobian Measurement model Jacobian with respect to the state
            %vector, x = [NED; Theta; vb; omega; eta].

            % Get necessary dimensions.
            ny = length(obj.Index);
            nx = size(x,1);
            nu = size(u,1);
            if isobject(Parameters)
                np = size(Parameters.Mean,1);
            else
                np = size(Parameters,1);
            end

            % Jacobians are constant.
            dhdx = zeros(ny,nx);
            dhdu = zeros(ny,nu);
            dhdp = zeros(ny,np);
            dhdx(1:ny,obj.Index) = eye(ny);
            
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
                newObj = StateEstimate(s.Name,s.Index); 
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