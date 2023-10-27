classdef LinearStaticStateFeedback < ControlLaw
    %LinearStaticStateFeedback Summary of this class goes here
    %   Detailed explanation goes here

    properties
        K
        x0
        u0
        FeedbackStateIndex
    end

    methods
        function obj = LinearStaticStateFeedback(K,x0,u0)
            %LinearStaticStateFeedback Construct an instance of this class
            %   Detailed explanation goes here
            obj.NumberOfStates = 0;
            obj.K = K;
            obj.x0 = x0;
            obj.u0 = u0;
        end

        function x_ctrldot = Dynamics(obj,~,~,~,~,~)
            x_ctrldot = zeros(0,1);
        end

        function u = Output(obj,~,~,x_rb,~,~)
            if isempty(obj.FeedbackStateIndex)
                x = x_rb;
            else
                x = x_rb(obj.FeedbackStateIndex,1);
            end
            dx = x - obj.x0;
            du = -obj.K*dx;
            u = du + obj.u0;
        end
    end
end