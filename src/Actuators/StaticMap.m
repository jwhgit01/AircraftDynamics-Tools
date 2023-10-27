classdef StaticMap < ActuatorModel
    %StaticMap Summary of this class goes here
    %   The polynomial should map raw digital units to physical units (i.e.
    %   PWM -> radians). The Bounds not only define the valid domain of
    %   this map, but also the trim point for which the physical units are
    %   zero. The independent variable of the polynomial is shifted by the
    %   trim value to cross through the origin. The inverse mapping
    %   polynomial is similarly shifted.

    properties
        Polynomial % up to 5th order (centered at trim value)
        InversePolynomial % linear approximation about trim
        Bounds % [min trim max]
    end

    methods

        function obj = StaticMap(InputAxis,Polynomial,Bounds)
            %StaticMap Construct an instance of this class
            obj.Type = 'StaticMap';
            obj.InputAxis = InputAxis;
            obj.Polynomial = Polynomial; 
            obj.Bounds = Bounds;

            % inverse mapping (of same order)
            % i.e. deg -> PWM
            n = length(Polynomial)-1;
            y = linspace(Bounds(1),Bounds(3),100).'-Bounds(2); % pwm points
            x = polyval(Polynomial,y); % PWM2deg at each point
            p = polyfit(x,y,n);
            obj.InversePolynomial = p;
        end % StaticMap

        function x_actdot = Dynamics(~,~)
            %Dynamics N/A for Static Map
            x_actdot = NaN;
        end

        function delta = Output(~,delta_cmd)
            %Output N/A for Static Map
            delta = delta_cmd;
        end

        function delta = Map(obj,u_cmd)
            %Map Convert an array in commanded units to actual units
            delta = zeros(size(u_cmd));
            for ii = 1:length(u_cmd)
                delta(ii) = polyval(obj.Polynomial,u_cmd(ii)-obj.Bounds(2));
            end
        end % Map

        function cmd = InverseMap(obj,delta)
            %Map Convert an array in actual units to commanded units
            cmd = zeros(size(delta));
            for ii = 1:length(delta)
                cmd(ii) = polyval(obj.InversePolynomial,delta(ii))+obj.Bounds(2);
            end
        end % InvMap
      
    end % public methods

    methods(Static)

        function obj = loadobj(s)
            %loadobj
            %
            % If any of the saved properties no longer are defined in the
            % class definition, create a dynamic property and give a
            % warning.
            if isstruct(s)
                newObj = StaticMap(s.InputAxis,s.Polynomial,s.Bounds); 
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
        % Import StaticMap factory   
            
            % CSV file
            CSVfile = ['AircraftData/' Aircraft.Name '_StaticMap.csv'];

            % check to see if file exists
            if 0 == exist(CSVfile,'file')
                error('Static Map CSV file does not exist.');
            end
            
            % import data as table
            tbl = readtable(CSVfile);

            % check that table is not empty
            if height(tbl) < 1
                error('Static Map CSV file is empty!');
            end
    
            % create StaticMap objects for each actuator
            for ii = 1:height(tbl)
                p = fliplr(tbl{ii,2:7}); % polynomial coefficients
                p(isnan(p)) = []; % remove empty orders
                ActuatorName = tbl{ii,1}{1}; % actuator name
                Bounds = tbl{ii,8:10}; % actuator bounds
                sm = StaticMap(ActuatorName,p,Bounds); % create object
                Aircraft.AddActuatorModel(sm,ii); % assign to aircraft
            end

        end % Import

    end % static methods


end % classdef