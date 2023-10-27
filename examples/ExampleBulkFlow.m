function W = ExampleBulkFlow(t,X)

    % Position components
    x = X(1,1);
    y = X(2,1);
    z = X(3,1);
    
    % Analytic Expression
%     Wx = 5*sin(y).*(1-sin(0.1*t));
%     Wy = cos(x).*atan(z+100);
%     Wz = -sin(0.1*x) - cos(0.1*y);

    Wx = 4;
    Wy = -3;
    Wz = -(1-cos(2*pi/200*x));

    % Wind vector
    W = [Wx; Wy; Wz];

end