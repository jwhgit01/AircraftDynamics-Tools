classdef (Abstract) Aircraft < handle & dynamicprops
%Aircraft Summary of this class goes here
%   Detailed explanation goes here

properties
    Name
    WorkspaceDir
    ToolboxVersion
    Mass
    MomentOfInertia
    ActuatorModel
    AerodynamicModel
    Sensors
    Units
end % properties

methods
    
    function obj = Aircraft(Name)
        %Aircraft Construct an empty instance of this class by name
        obj.Name = Name;
        
        % save the full path of this aircraft's workspace
        obj.WorkspaceDir = [pwd filesep];
    end

    function Save(obj,AircraftPath)
        %Save Save an aircraft model.
        if nargin < 2
            filepath = obj.WorkspaceDir;
            filename = [filepath filesep obj.Name '.mat'];
        else
            filename = AircraftPath;
        end
        evalin('base',['save(''' filename ''',''' obj.Name ''',''-v7.3'');']);   
    end % Save

    function b = saveobj(a)
        %saveobj
        %
        % When saving an aircraft object, set the version.
        AircraftClassDir = which('Aircraft.m');
        [filepath,~,~] = fileparts(AircraftClassDir);
        VersionDir = [filepath filesep '..' filesep '..' filesep 'VERSION'];
        ver = fileread(VersionDir);
        a.ToolboxVersion = ver;
        b = a;
    end

    function ImportInertialParameters(obj)
        %SetAircraftParameters set the aircraft parameters
        
        % check to see if file exists
        if 0 == exist(['AircraftData/' obj.Name '_InertialParameters.csv'],'file')
            error('InertialParameters CSV file does not exist.');
        end
        
        % read csv data as table
        tbl = readtable(['AircraftData/' obj.Name '_InertialParameters.csv']);
        
        % check that data is there
        if height(tbl) < 1
            error('InertialParameters.csv empty!');
        end

        % assign to aircraft object
        obj.Mass = tbl.mass_kg(1);
        Ixx = tbl.Ixx_kgm2(1);
        Iyy = tbl.Iyy_kgm2(1);
        Izz = tbl.Izz_kgm2(1);
        Ixz = tbl.Ixz_kgm2(1);
        obj.MomentOfInertia = [Ixx,0,-Ixz;0,Iyy,0;-Ixz,0,Izz];
    
    end % ImportInertialParameters

    
    function AddActuatorModel(obj,ActModelObj,InputNumber)
        %AddActuatorModel add an actuator model to the aircraft object
        obj.ActuatorModel.(ActModelObj.Type)(InputNumber,1) = ActModelObj;
    end % AddActuatorModel

    function AddAerodynamicModel(obj,AeroModelObj,Name)
        %AddAerodynamicModel add an aerodynamic model to the aircraft object
        obj.AerodynamicModel.(Name) = AeroModelObj;
    end % AddAerodynamicModel

    function AddSensor(obj,SensorObj)
        %AddSensor add a sensor to the aircraft object
        mc = metaclass(SensorObj);
        obj.Sensors.(mc.Name).(SensorObj.Name) = SensorObj;
    end % AddSensor
    
end % public methods

methods(Static)

    function CreateWorkspace(AircraftName)
        % CreateWorkspace Create a new aircraft workspace and object
        
        % check to see if this aircraft already exists
        if 7 == exist(AircraftName,'dir')
            error('Aircraft workspace already exists!');
        end
        
        % create the folder structure for a system identification workspace
        mkdir(AircraftName);
        mkdir([AircraftName '/AircraftData']);
        
        % copy template files
        temp = which('AircraftDynamcis-Tools/lib/InertialParameters.csv');
        [libpath,~,~] = fileparts(temp);
        copyfile([libpath filesep 'InertialParameters.csv'],[AircraftName '/AircraftData/' AircraftName '_InertialParameters.csv']);
        copyfile([libpath filesep 'StaticMap.csv'],[AircraftName '/AircraftData/' AircraftName '_StaticMap.csv']);
        copyfile([libpath filesep 'TimeDelay.csv'],[AircraftName '/AircraftData/' AircraftName '_TimeDelay.csv']);
        copyfile([libpath filesep 'RateLimit.csv'],[AircraftName '/AircraftData/' AircraftName '_RateLimit.csv']);
        copyfile([libpath filesep 'FirstOrder.csv'],[AircraftName '/AircraftData/' AircraftName '_FirstOrder.csv']);
        copyfile([libpath filesep 'AerodynamicParameters.csv'],[AircraftName '/AircraftData/' AircraftName '_AerodynamicParameters_NAME.csv']);
        
    end % CreateWorkspace

    function Load(AircraftPath)
        %Load Load an aircraft model into the workspace and create a
        % backup just in case.
        
        % copy the .mat file to a backup
        [filepath,AircraftName,~] = fileparts(AircraftPath);
        copyfile(AircraftPath,[filepath filesep AircraftName '_bk.mat']);
        
        % load the aircraft into the workspace
        load(AircraftPath,AircraftName);
        eval(['assignin(''base'',AircraftName,' AircraftName ');']);
        
    end % Load

    function [y,H] = MeasurementModel(x,u,Parameters,Constants,ideriv,ListOfSensors,ParameterEstimation)
        %MeasurementModel Compute the modeled outputs of the aircraft's
        % sensors listed in ListOfSensors. If ListOfSensors equals
        % 'FullState', then the sensor model is the state vector,
        % x = [NED; Theta; vb; omega; eta].
        
        % check if full state output or partial state
        if strcmp(ListOfSensors,'FullState')
            y = x;
            H = eye(size(x,1));
            return
        end
        if isnumeric(ListOfSensors)
            y = x(ListOfSensors,1);
            H = zeros(length(ListOfSensors),size(x,1));
            H(ListOfSensors,ListOfSensors) = eye(length(ListOfSensors));
            return
        end
        
        nz = Constants.nz;
        ns = length(ListOfSensors);
        nx = size(x,1);
        y = zeros(nz,1);
        if ideriv
            if ParameterEstimation
                H = zeros(nz,nx+length(Parameters));
            else
                H = zeros(nz,nx);
            end
        else
            H = [];
        end
        idx0 = 0;
        for ii = 1:ns
            Sen = ListOfSensors{ii};
            yi = Sen.Model(x,u,Parameters,Constants);
            idx = (1:size(yi,1)) + idx0;
            y(idx,1) = yi;
            if ideriv
                [dhdx,~,dhdp] = Sen.Jacobian(x,u,Parameters,Constants);
                if ParameterEstimation
                    H(idx,:) = [dhdx,dhdp];
                else
                    H(idx,:) = dhdx;
                end
            end
            idx0 = idx(end); 
        end

    end % MeasurementModel

end % static methods

methods(Abstract)
    obj = loadobj(s)
end % abstract methods

methods(Abstract,Static)
        NewAircraft(AircraftName)
    end

end % classdef