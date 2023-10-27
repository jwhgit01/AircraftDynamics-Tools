close all
clear
clc
LaTeXify

% Create the simulation enviroment
Name = 'MTD2_Sim_WingsLevel_UniformPlusDryden_v1'
Simulation = AircraftSimulator(Name);

% Set the sim time
Simulation.FinalTime = 30; 
N = length(Simulation.Time);

% Load an aircraft
Aircraft.Load('G:\My Drive\Research\AircraftModels\MTD2\MTD2.mat')
Simulation.Aircraft = MTD2;

% Create an envirment model
Enviroment = EnviromentModel;
Enviroment.BulkFlowField = @(t,X) [1;-3;-2];
Enviroment.BulkGradient = []; % Ignore bulk gradient effects in simulation

% Create a Dryden turbulence model
Altitude = 100;
Specification = 'MIL-HDBK-1797B';
Windspeed = 10;
Airspeed = 18;
Span = MTD2.Span;
N_nu = floor(N/10);
nu = randn(N_nu,1); % band-limited
t_nu = linspace(0,Simulation.FinalTime,N_nu).';
WhiteNoiseFunction = @(t) interp1(t_nu,nu,t,'previous',0);
Dryden = DrydenTurbulence(Altitude,Specification,Windspeed,Airspeed,Span,WhiteNoiseFunction);
Dryden.Simulation = Simulation;
Enviroment.TurbulenceModel = Dryden;

% Set the air density and gravity to be constant
[~,~,~,rho] = atmosisa(500);
Enviroment.AirDensity = rho;
Enviroment.Gravity = gravitywgs84(500,37.19711);

% Add enviroment model to simulation
Simulation.EnviromentModel = Enviroment;

% Select aerodynamic model and parameter estimates
AeroModel = MTD2.AerodynamicModel.WingsLevelCruisePlusThrust;
Simulation.AerodynamicModel = AeroModel;
AeroParams = Simulation.AerodynamicModel.ParameterEstimates.WingsLevelCruisePlusThrust;
Simulation.AerodynamicParameters = AeroParams;

% Use Quaternions as the attitude parameterization
Simulation.AttitudeParameterization = 'EulerAngles';

% Compute the linearization of the system about arbitrary trim
NED0 = zeros(3,1);
alpha0 = 1;
Theta0 = [0; alpha0*pi/180; 0];
vb0 = [18*cosd(alpha0);0;-18*sind(alpha0)];
omega0 = zeros(3,1);
delta0 = [0; 0; 0; 220];
x_rb_0 = [NED0;Theta0;vb0;omega0];
I = Simulation.Aircraft.MomentOfInertia;
g = Enviroment.Gravity;
Constants.rho = rho;
A_rb = RigidBodyJacobian(x_rb_0,I,g,Simulation.AttitudeParameterization);
[dfdv,dfdomega,dfddelta,~] = AeroModel.Jacobian(zeros(2,1),vb0,omega0,delta0,AeroParams,Constants);
A_aero = [zeros(6,12);zeros(6,6),dfdv,dfdomega];
A12 = A_rb + A_aero;
B12 = [zeros(6,4); dfddelta];

% Create an open loop control law
% ControlFunction = @(t) [0;-1*pi/180;0;220];
% Controller = OpenLoop(ControlFunction);
% Simulation.ControlLaw = Controller;

% % Create a simple LQR control law
FeedbackStateIndex = 4:12;
xnom = x_rb_0(FeedbackStateIndex,1);
A = A12(FeedbackStateIndex,FeedbackStateIndex);
B = B12(FeedbackStateIndex,:);
dTheta_max = 5*pi/180;
dvb_max = 0.5;
domega_max = 30*pi/180;
ddelta_max = 10*pi/180;
dOmega_max = 30;
Q = 10*diag([1/(dTheta_max^2)*ones(1,3), 1/(dvb_max^2)*ones(1,3), 1/(domega_max^2)*ones(1,3)]);
R = diag([1/(ddelta_max^2)*ones(1,3), 1/(dOmega_max^2)]);
K = lqr(A,B,Q,R);
Controller = LinearStaticStateFeedback(K,xnom,delta0);
Controller.FeedbackStateIndex = FeedbackStateIndex;
Simulation.ControlLaw = Controller;

% Simple H2 controller
% dTheta_max = 5*pi/180;
% dvb_max = 1;
% domega_max = 30*pi/180;
% ddelta_max = 10*pi/180;
% dOmega_max = 30;
% w_max = 5;
% FeedbackStateIndex = 4:12;
% xnom = x_rb_0(FeedbackStateIndex,1);
% A = A12(FeedbackStateIndex,FeedbackStateIndex);
% Bu = B12(FeedbackStateIndex,:);
% Bw = w_max*A_aero(FeedbackStateIndex,7:9);
% Q = diag([1/(dTheta_max^2)*ones(1,3), 1/(dvb_max^2)*ones(1,3), 1/(domega_max^2)*ones(1,3)]);
% R = 0.001*diag([ddelta_max^2*ones(1,3), dOmega_max^2]);
% Cy = eye(9);
% C = [zeros(9,3), Cy];
% Dyu = zeros(9,4);
% Dyw = zeros(9,3);
% Cz = [chol(Q); zeros(4,9)];
% Dzu = [zeros(9,4); chol(R)];
% Dzw = zeros(9+4,3);
% Ba = [Bw Bu];
% Ca = [Cz; Cy];
% Da = [Dzw, Dzu; Dyw, Dyu];
% P = ss(A,Ba,Ca,Da);
% [K,~,gamma] = h2syn(P,9,4)
% Controller = LinearDynamicOutputFeedback(C,K,x_rb_0,delta0);
% Simulation.ControlLaw = Controller;

% Set the initial condition
x0 = InitialCondition(Simulation);
x0.RigidBody.Position(3,1) = -Altitude;
x0.RigidBody.Attitude = [0;0;0];
x0.RigidBody.Velocity = [18*cosd(3);0;-18*sind(3)];
%
w0 = Enviroment.BulkFlowField(0,x0.RigidBody.Position);
vr0 = x0.RigidBody.Velocity - w0;
alphabeta0 = [atan2(vr0(3,1),vr0(1,1));asin(vr0(2,1)/norm(vr0))];
x0.Unsteady = AeroModel.DerivativeTimeConstant*alphabeta0;

% Run the simulation
opts = odeset('RelTol',1e-3,'AbsTol',1e-6,'OutputFcn',@odeplot);
[t,x] = ode23t(@Simulation.Dynamics,Simulation.Time,x0.IC,opts);

% Parse the state vector time history
[x_rb,x_us,x_act,x_env,x_obsv,x_ctrl] = Simulation.State(x);

% Compute the input and wind time histories
u = zeros(N,4);
w = zeros(N,3);
for k = 1:N
    x_ctrlk = x_ctrl(k,:).';
    x_rbk = x_rb(k,:).';
    x_envk = x_env(k,:).';
    x_obsvk = x_obsv(k,:);
    u(k,:) = Controller.Output(t(k),x_ctrlk,x_rbk,x_obsvk,[]).';
    [wk,~] = Enviroment.Output(t(k),x_envk,x_rbk);
    w(k,:) = wk.';
end

% Plot
N = x(:,1);
E = x(:,2);
D = x(:,3);
if strcmp(Simulation.AttitudeParameterization,'Quaternion')
    EulerAngles = fliplr(quat2eul(x(:,4:7),'ZYX'));
elseif strcmp(Simulation.AttitudeParameterization,'EulerAngles')
    EulerAngles = x(:,4:6);
else
    error('TODO')
end
phi = EulerAngles(:,1);
theta = EulerAngles(:,2);
psi = EulerAngles(:,3);
scaleFactor = 10;
NumAircraft = 20;
aircraft = 'cessna';
figure(1)
clf
trajectory(N,E,D,phi,theta,psi,scaleFactor,NumAircraft,aircraft)

% plot inputs
figure(2)
plot(t,u(:,1:3)*180/pi)
ylim([-30 30])
ylabel('Actuator deflection [deg]')
yyaxis right
plot(t,u(:,4))
ylim([100 300])
ylabel('Propeller Speed [rad/s]')
legend('$\delta a$','$\delta e$','$\delta r$','$\Omega$')

% plot wind
figure(3)
plot(t,w)
ylabel('Wind Velocity [m/s]')
legend('$w_x$','$w_y$','$w_z$')


%% Save Restuls
vars = {'Simulation','thist','x','u','w','x_rb','x_us','x_act','x_env','x_obsv','x_ctrl'};
matfile = Name
