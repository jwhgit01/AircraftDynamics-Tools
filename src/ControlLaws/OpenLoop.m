classdef OpenLoop < ControlLaw
    %OpenLoop Summary of this class goes here
    %   Detailed explanation goes here

    properties
        ControlFunction
    end

    methods
        function obj = OpenLoop(ControlFunction)
            %OpenLoop Construct an instance of this class
            %   Detailed explanation goes here
            obj.ControlFunction = ControlFunction;
            obj.NumberOfStates = 0;
        end

        function x_ctrldot = Dynamics(obj,~,~,~,~,~)
            x_ctrldot = zeros(0,1);
        end

        function u = Output(obj,t,~,~,~,~)
            u = obj.ControlFunction(t);
        end
    end
end