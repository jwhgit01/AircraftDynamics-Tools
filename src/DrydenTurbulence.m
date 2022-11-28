function H = DrydenTurbulence(Altitude,Specification,WindSpeed,Airspeed,Span)
% 
% Source: https://www.mathworks.com/help/aeroblks/drydenwindturbulencemodelcontinuous.html

% Determine characteristic lengths, L for low altitude turbulence (<1000ft)
if strcmp(Specification,'MIL-F-8785C')
    Lu = Altitude/((0.177+0.000823*Altitude)^1.2);
    Lv = Lu;
    Lw = Altitude;
elseif strcmp(Specification,'MIL-HDBK-1797B')
    Lu = Altitude/((0.177+0.000823*Altitude)^1.2);
    Lv = Lu/2;
    Lw = Altitude/2;
else
    error('Invalid Specification');
end

% Turbulence standard deviations using low altitude model
sigma_w = 0.1*WindSpeed;
sigma_u = sigma_w/((0.177+0.000823*Altitude)^1.2);
sigma_v = sigma_u;

% continuous time transfer functions
if strcmp(Specification,'MIL-F-8785C')
    
    % Longitudinal
    H.u = tf(sigma_u*sqrt(2*Lu/(pi*Airspeed)),[Lu/Airspeed,1]);
    H.p = tf(sigma_w*sqrt(0.8/Airspeed)*(pi/(4*Span))^(1/6),Lw^(1/3)*[4*Span/(pi*Airspeed),1]);

    % Lateral
    H.v = tf(sigma_v*sqrt(Lv/(pi*Airspeed))*[sqrt(3)*Lv/Airspeed,1],[(Lv/Airspeed)^2,2*Lv/Airspeed,1]);
    H.r = tf([1/Airspeed,0],[3*Span/(pi*Airspeed),1])*H.v;

    % Vertical
    H.w = tf(sigma_w*sqrt(Lw/(pi*Airspeed))*[sqrt(3)*Lw/Airspeed,1],[(Lw/Airspeed)^2,2*Lw/Airspeed,1]);
    H.q = tf([1/Airspeed,0],[4*Span/(pi*Airspeed),1])*H.w;

elseif strcmp(Specification,'MIL-HDBK-1797B')
    
    % Longitudinal
    H.u = tf(sigma_u*sqrt(2*Lu/(pi*Airspeed)),[Lu/Airspeed,1]); % same as MIL-F-8785C
    H.p = tf(sigma_w*sqrt(0.8/Airspeed)*(pi/(4*Span))^(1/6),(2*Lw)^(1/3)*[4*Span/(pi*Airspeed),1]);

    % Lateral
    H.v = tf(sigma_v*sqrt(2*Lv/(pi*Airspeed))*[2*sqrt(3)*Lv/Airspeed,1],[(2*Lv/Airspeed)^2,4*Lv/Airspeed,1]);
    H.r = tf([1/Airspeed,0],[3*Span/(pi*Airspeed),1])*H.v; % same form as MIL-F-8785C

    % Vertical
    H.w = tf(sigma_w*sqrt(2*Lw/(pi*Airspeed))*[2*sqrt(3)*Lw/Airspeed,1],[(2*Lw/Airspeed)^2,4*Lw/Airspeed,1]);
    H.q = tf([1/Airspeed,0],[4*Span/(pi*Airspeed),1])*H.w; % same form as MIL-F-8785C

else
    error('Invalid Specification');
end

end