classdef EnviromentModel < handle
%EnviromentModel 
%
% Copyright (c) 2023 Jeremy W. Hopwood. All rights reserved.
%
% This class ...
%

properties
    BulkFlowField % constant or function of position or time and position
    BulkGradient % constant or function of position or time and position
    TurbulenceModel % VonKarman, Dryden, Kaimal, GaussMarkov
    AirDensity % Constant or 'ISA'
    Gravity % Constant or 'WGS84'
end
    
methods

    function obj = EnviromentModel()
        %EnviromentModel Create an instance with default properties
        obj.BulkFlowField = @(t,x) zeros(3,1);
        obj.BulkGradient = @(t,x) zeros(3,3);
        obj.AirDensity = 1.225;
        obj.Gravity = 9.81;
    end

    function [w,gradW] = Output(obj,t,x_env,x_rb)
        %Output Compute the wind velocity and gradient

        % Compute the output of the turbulence model if applicable
        if isempty(obj.TurbulenceModel)
            dw = zeros(3,1);
            dgradW = zeros(3,3);
        else
            [dw,dgradW] = obj.TurbulenceModel.Output(t,x_env,x_rb);
        end

        % Evaluate bulk flow
        NED = x_rb(1:3,1);
        if isempty(obj.BulkFlowField)
            wbar = zeros(3,1);
        else
            wbar = obj.BulkFlowField(t,NED);
        end
        if isempty(obj.BulkGradient)
            gradWbar = zeros(3,3);
        else
            gradWbar = obj.BulkGradient(t,NED);
        end

        % Add the turbulence to the bulk flow
        w = wbar + dw;
        gradW = gradWbar + dgradW;

    end

end % public methods

end