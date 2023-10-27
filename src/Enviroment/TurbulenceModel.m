classdef (Abstract) TurbulenceModel < handle
%DrydenTurbulence 
%
% Copyright (c) 2023 Jeremy W. Hopwood. All rights reserved.
%
% This class ...
%


properties
    NumberOfStates
end
    
methods (Abstract)
    x_envdot = Dynamics(obj,t,x_env,x_rb)
    [dw,dgradW] = Output(obj,t,x_env,x_rb)
end

end