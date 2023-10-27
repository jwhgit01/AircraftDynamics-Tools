classdef FixedWingAircraft < Aircraft
    %FixedWing Summary of this class goes here
    %   Detailed explanation goes here

    properties
        Span
        Chord
        WingArea
        PropDiameter
    end % properties

    methods

        function obj = FixedWingAircraft(Name)
            %FixedWing Construct an instance of this class
            %   Call Aircraft constructor
            obj@Aircraft(Name); 
        end

        function ImportFixedWingGeometry(obj)
            %ImportFixedWingGeometry set the aircraft parameters
            
            % check to see if file exists
            if 0 == exist(['AircraftData/' obj.Name '_FixedWingGeometry.csv'],'file')
                error('InertialParameters CSV file does not exist.');
            end
            
            % read csv data as table
            tbl = readtable(['AircraftData/' obj.Name '_FixedWingGeometry.csv']);
            
            % check that data is there
            if height(tbl) < 1
                error('FixedWingGeometry CSV empty!');
            end

            % assign to aircraft object
            obj.Span = tbl.span_m(1);
            obj.Chord = tbl.chord_m(1);
            obj.WingArea = tbl.wingarea_m2(1);

        end % ImportFixedWingGeometry        

    end % public methods

    methods(Static)

        function obj = loadobj(s)
            %loadobj
            %
            % If any of the saved properties no longer are defined in the
            % class definition, create a dynamic property and give a
            % warning.
            if isstruct(s)
                newObj = FixedWingAircraft(s.Name); 
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
            % NewAircraft Create a new aircraft workspace and FixedWing object
            
            % General Initialization
            Aircraft.CreateWorkspace(AircraftName);

            % FixedWing-specific initialization
            temp = which('AircraftDynamcis-Tools/lib/InertialParameters.csv');
            [libpath,~,~] = fileparts(temp);
            copyfile([libpath filesep 'FixedWingGeometry.csv'],[AircraftName '/AircraftData/' AircraftName '_FixedWingGeometry.csv']);

            % create FixedWingAircraft object
            eval([AircraftName '=FixedWingAircraft(AircraftName);']);
            
            % save object to mat file
            save([AircraftName filesep AircraftName '.mat'], AircraftName);
            
            % sucess!
            disp([AircraftName ' created successfully!']);

        end % NewAircraft

    end % static methods

end % classdef