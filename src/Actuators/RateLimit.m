classdef RateLimit < ActuatorModel
    %RateLimit Summary of this class goes here
    %   Detailed explanation goes here

    properties
        RateLim % note: can be Inf
    end

    methods

        function obj = RateLimit(InputAxis,RateLim)
            %StaticMap Construct an instance of this class
            obj.Type = 'RateLimit';
            obj.InputAxis = InputAxis;
            obj.RateLim = RateLim;
        end

        function deltadot = Dynamics(delta,delta_cmd)
            %Dynamics Define the acutator dynamics
            deltadot = -sat(delta-delta_cmd)*obj.RateLim;
        end

        function u = Map(obj,t,cmd)
            %Map Convert an array of commanded values to actual values
            u = cmd;
            for ii = 2:length(u)
                dudt = (u(ii)-u(ii-1))/(t(ii)-t(ii-1));
                if abs(dudt) > obj.RateLim
                    u(ii) = u(ii-1) + sign(dudt)*obj.RateLim*(t(ii)-t(ii-1));
                end
            end
        end

        % note the inverse mapping is non-unique and thus N/A
      
    end % public methods

    methods (Access = private)
        function y = sat(x)
            xoverepsil = x*1e6;
            if abs(xoverepsiln) <= 1
                y = xoverepsil;
            else
                y = sign(xoverepsil);
            end
        end
    end % private methods

    methods(Static)

        function obj = loadobj(s)
            %loadobj
            %
            % If any of the saved properties no longer are defined in the
            % class definition, create a dynamic property and give a
            % warning.
            if isstruct(s)
                newObj = RateLimit(s.InputAxis,s.RateLim); 
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

        function Import(Aircraft)
        % Import RateLimit factory   

            % CSV file
            CSVfile = ['AircraftData/' Aircraft.Name '_RateLimit.csv'];

            % check to see if file exists
            if 0 == exist(CSVfile,'file')
                error('Rate Limit CSV file does not exist.');
            end
            
            % import data as table
            tbl = readtable(CSVfile);

            % check that table is not empty
            if height(tbl) < 1
                error('Rate Limit CSV file is empty!');
            end
    
            % create RateLimit objects for each actuator
            for ii = 1:height(tbl)
                rl_rad_s = tbl{ii,2}; % rate limit [rad/s]
                ActuatorName = tbl{ii,1}{1}; % actuator name
                rl = RateLimit(ActuatorName,rl_rad_s); % create object
                Aircraft.AddActuatorModel(rl); % assign to aircraft
            end

        end % Import

    end % static methods

end % classdef