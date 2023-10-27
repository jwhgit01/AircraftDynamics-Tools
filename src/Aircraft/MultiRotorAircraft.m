classdef MultiRotorAircraft < Aircraft
    %MultiRotorAircraft Summary of this class goes here
    %   Detailed explanation goes here

    properties
        RotorRadius
        NumberOfRotors
        ArmLength
        HubHeight
        MixingMatrix
        MotorMomentOfInertia % scalar (about z-axis)
    end % properties

    methods

        function obj = MultiRotorAircraft(Name)
            %MultiRotorAircraft Construct an instance of this class
            %   Call Aircraft constructor
            obj@Aircraft(Name); 
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
                newObj = MultiRotorAircraft(s.Name); 
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

        function NewAircraft(AircraftName)
            % NewAircraft Create a new aircraft workspace and MultiRotorAircraft object
            
            % General Initialization
            Aircraft.CreateWorkspace(AircraftName);

            % create MultiRotorAircraft object
            eval([AircraftName '=MultiRotorAircraft(AircraftName);']);
            
            % save object to mat file
            save([AircraftName filesep AircraftName '.mat'], AircraftName);
            
            % sucess!
            disp([AircraftName ' created successfully!']);

        end % NewAircraft

    end % static methods

end % classdef