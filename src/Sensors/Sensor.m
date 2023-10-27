classdef (Abstract) Sensor < handle & dynamicprops
    %Sensor Abstract class for actuator models
    %   Example subclasses:
    %       Accelerometer
    %       Gyroscope
    %       Magnetometer
    %       AirDataUnit
    %       PitotTube
    %       GPS

    properties
        Name
        Covariance
        ULOGTopic
    end % properties

    methods(Abstract)
        y = Model(obj,x,u,Parameters,Constants)
        [dhdx,dhdu,dhdp] = Jacobian(obj,x,u,Parameters,Constants)
        obj = loadobj
    end % Abstract methods

end % classdef