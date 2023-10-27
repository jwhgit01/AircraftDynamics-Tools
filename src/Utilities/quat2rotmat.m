function R = quat2rotmat(q)
%quat2rotmat
%
% Copyright (c) 2023 Jeremy W. Hopwood. All rights reserved.
%
% This function ...
%
% Inputs:
%
%   in1     The first input...
%  
% Outputs:
%
%   out1    The first output...
%

qr = q(1,1);
qi = q(2,1);
qj = q(3,1);
qk = q(4,1);
s = norm(q);

R = [1-2*s*(qj^2+qk^2),2*s*(qi*qj-qk*qr),2*s*(qi*qk+qj*qr);
     2*s*(qi*qj+qk*qr),1-2*s*(qi^2+qk^2),2*s*(qj*qk-qi*qr);
     2*s*(qi*qk-qj*qr),2*s*(qj*qk+qi*qr),1-2*s*(qi^2+qj^2)];

end