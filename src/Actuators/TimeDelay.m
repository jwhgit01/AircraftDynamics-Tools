classdef TimeDelay < ActuatorModel
    %TimeDelay Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Delay % [s]
    end

    methods
        function obj = TimeDelay(InputAxis,Delay)
            %TimeDelay Construct an instance of this class
            obj.Type = 'TimeDelay';
            obj.InputAxis = InputAxis;
            obj.Delay = Delay; 
        end

        function deltadot = Dynamics(~,~)
            %Dynamics N/A for Time Delay
            deltadot = NaN;
        end

        function u_out = Map(obj,t_cmd,u_cmd)
            %Map shift the commanded value by the time delay
            t_shift = t_cmd + obj.Delay;
            u_out = zeros(size(u_cmd));
            for ii = 1:length(t_cmd)
                u_out(ii) = interp1(t_shift,u_cmd,t_cmd(ii));
            end
        end

        function u_cmd = InvMap(obj,t_out,u_out)
            %Map shift the delayed value back to the commanded
            t_shift = t_out - obj.Delay;
            u_cmd = zeros(size(u_out));
            for ii = 1:length(t_out)
                u_cmd(ii) = interp1(t_shift,u_out,t_cmd(ii));
            end
        end
      
    end % public methods

    methods(Static)

        function obj = loadobj(s)
            %loadobj
            %
            % If any of the saved properties no longer are defined in the
            % class definition, create a dynamic property and give a
            % warning.
            if isstruct(s)
                newObj = TimeDelay(s.InputAxis,s.Delay); 
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
        % Import TimeDelay factory   

            % CSV file
            CSVfile = ['AircraftData/' Aircraft.Name '_TimeDelay.csv'];

            % check to see if file exists
            if 0 == exist(CSVfile,'file')
                error('Time Delay CSV file does not exist.');
            end
            
            % import data as table
            tbl = readtable(CSVfile);

            % check that table is not empty
            if height(tbl) < 1
                error('Time Delay CSV file is empty!');
            end
    
            % create TimeDelay objects for each actuator
            for ii = 1:height(tbl)
                tau = tbl{ii,2}; % delay [s]
                ActuatorName = tbl{ii,1}{1}; % actuator name
                td = TimeDelay(ActuatorName,tau); % create object
                Aircraft.AddActuatorModel(td); % assign to aircraft
            end

        end % Import

    end % static methods

end % classdef