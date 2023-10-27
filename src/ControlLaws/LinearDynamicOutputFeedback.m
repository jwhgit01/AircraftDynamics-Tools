classdef LinearDynamicOutputFeedback < ControlLaw
    %LinearDynamicOutputFeedback Summary of this class goes here
    %   Detailed explanation goes here

    properties
        C
        K % ss object
        x0
        u0
    end

    methods
        function obj = LinearDynamicOutputFeedback(C,K,x0,u0)
            %LinearDynamicOutputFeedback Construct an instance of this class
            %   Detailed explanation goes here
            obj.NumberOfStates = size(K.A,1);
            obj.C = C;
            obj.K = K;
            obj.x0 = x0;
            obj.u0 = u0;
        end

        function x_ctrldot = Dynamics(obj,~,x_ctrl,x_rb,~,~)
            dx = x_rb - obj.x0;
            dy = obj.C*dx;
            x_ctrldot = obj.K.A*x_ctrl + obj.K.B*dy;
        end

        function u = Output(obj,~,x_ctrl,x_rb,~,~)
            dx = x_rb - obj.x0;
            dy = obj.C*dx;
            du = obj.K.C*x_ctrl + obj.K.D*dy;
            u = du + obj.u0;
        end
    end
end