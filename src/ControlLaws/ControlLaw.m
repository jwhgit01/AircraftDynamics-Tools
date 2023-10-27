classdef (Abstract) ControlLaw < handle
    %ControlLaw Summary of this class goes here
    %   Detailed explanation goes here

    properties
        NumberOfStates
        Parameters
        Reference
    end

    methods (Abstract)
        x_ctrldot = Dynamics(obj,t,x_ctrl,x_rb,x_obsv,y);
        u = Output(obj,t,x_ctrl,x_rb,x_obsv,y);
    end
end