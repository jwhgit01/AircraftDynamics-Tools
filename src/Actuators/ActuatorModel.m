classdef (Abstract) ActuatorModel < handle & dynamicprops
    %ActuatorModel Abstract class for actuator models
    %   Example subclasses:
    %       StaticMap
    %       TimeDelay
    %       RateLimit
    %       FirstOrder

    properties
        Type % string that equals the subclass name
        InputAxis % string
        NumberOfStates
    end
    
    methods(Abstract)

        % Define the acutator dynamics if applicable. If not set to NaN
        x_actdot = Dynamics(x_act,delta_cmd)

        % Define the acutator output if applicable. If not set to NaN
        delta = Output(x_act,delta_cmd)
        
        % Map time series of known input values to true actuator values
        u_out = Map(obj,t_cmd,u_cmd)

        % If any of the saved properties no longer are defined in the
        % class definition, create a dynamic property and give a warning.
        obj = loadobj(s)

    end

    methods(Abstract = true, Static = true)
        Import(Aircraft) % import from CSV and assign to Aircraft
    end

end % classdef