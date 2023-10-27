classdef FirstOrder < ActuatorModel
    %FirstOrder Summary of this class goes here
    %   Detailed explanation goes here

    properties
        TimeConstant % [s]
    end

    methods
        function obj = FirstOrder(InputAxis,TimeConstant)
            %FirstOrder Construct an instance of this class
            obj.Type = 'FirstOrder';
            obj.InputAxis = InputAxis;
            obj.TimeConstant = TimeConstant; 
            obj.NumberOfStates = 1;
        end

        function x_actdot = Dynamics(x_act,delta_cmd)
            %deltadot Define the acutator dynamics
            tau = obj.TimeConstant;
            x_actdot = -(1/tau)*(x_act-delta_cmd);
        end

        function delta = Output(x_act,~)
            %deltadot Define the acutator output
            delta = x_act;
        end

        function u_out = Map(obj,t_cmd,u_cmd)
            %Map apply first-order map
            sys = tf(1,[obj.TimeConstant 1]);
            u_out = lsim(sys,u_cmd,t_cmd,u_cmd(1,:));
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
                newObj = FirstOrder(s.InputAxis,s.TimeConstant); 
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
        % Import FirstOrder factory   

            % CSV file
            CSVfile = ['AircraftData/' Aircraft.Name '_FirstOrder.csv'];

            % check to see if file exists
            if 0 == exist(CSVfile,'file')
                error('First Order CSV file does not exist.');
            end
            
            % import data as table
            tbl = readtable(CSVfile);

            % check that table is not empty
            if height(tbl) < 1
                error('First Order CSV file is empty!');
            end
    
            % create FirstOrder objects for each actuator
            for ii = 1:height(tbl)
                tau = tbl{ii,2}; % time constant [s]
                ActuatorName = tbl{ii,1}{1}; % actuator name
                fo = FirstOrder(ActuatorName,tau); % create object
                Aircraft.AddActuatorModel(fo,ii); % assign to aircraft
            end

        end % Import

    end % static methods

end % classdef