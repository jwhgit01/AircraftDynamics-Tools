function [Xequil,Uequil,ay,az] = TrimmingAlgorithm(model,Xguess,Uguess,orient,cost,acft)
%     1. Wings Level (gamma = 0)
%     2. Wings Level (gamma <> 0)
%     3. Steady Constant Altitude Turn
%     4. Steady Pull Up
%     5. Steady climbing/descending turn
%  The user will also be prompted for the number of iterations to be used in the numerical
%  minimization search.

format long

if(nargin==5)
xg=Xguess;
ug=Uguess;
else
xg=zeros(12,1);
ug=zeros(4,1);
end
%        gamma singam rr  pr   tr  phi cphi sphi thetadot coord stab  orient
const = [0.0    0.0   0.0 0.0  0.0 0.0 1.0  0.0   0.0      0.0   0.0  1];
rtod = 57.295779513082323;

% orient = menu('Choose an A/C Orientation','Wings Level (gamma = 0)',...
% 'Wings Level (gamma <> 0)','Steady Turn','Steady Pull Up');
const(12) = orient;
ndof = 6;

if orient == 1
   xg(1) = input('Velocity Vector (VT): ');
   xg(12) = input('Altitude (h): ');
end

if orient == 2
   xg(1) = 1.68781*input('Velocity Vector (kts) (VT): ');
   xg(12) = input('Altitude (ft) (h): ');
   gamm = input('Gamma (deg): ');
   const(1) = gamm/rtod;
   const(2) = sin(const(1));
end

if orient == 3
   xg(1) = 1.68781*input('Velocity Vector (kts) (VT): ');
   xg(12) = input('Altitude (ft) (h): ');
   psidot = input('Turn Rate (deg/s) (Psi dot): ');
   const(5) = psidot/rtod;
end

if orient == 4
   xg(1) = 1.68781*input('Velocity Vector (kts) (VT): ');
   xg(12) = input('Altitude (ft) (h): ');
   thetadot = input('Pitch Rate (deg/s) (Theta dot): ');
   const(9) = thetadot/rtod;
end

if orient == 5
   xg(1) = 1.68781*input('Velocity Vector (kts) (VT): ');
   xg(12) = input('Altitude (ft) (h): ');
   psidot = input('Turn Rate (deg/s) (Psi dot): ');
   const(5) = psidot/rtod;
   gamm = input('Gamma (deg): ');
   const(1) = gamm/rtod;
   const(2) = sin(const(1));
end

% Set up the initial guess for the state and control vectors

% if nargin~=5
% disp(' ')
% disp('Next Input The Initial Guess For The Equilibrium State And Control Vectors')
% disp('Remember To Match The Altitude and Air Speed You Just Keyed In:')
% disp(' ')
% getinput
% end

yesno = 1;
clear s

if orient == 3 %steady constant alt turn
s(1)=ug(1);
s(2)=ug(2);
s(3)=ug(3);
s(4)=ug(4);
s(5)=xg(2);
s(6)=xg(4);
s(7)=xg(5);
else
s(1)=ug(1);
s(2)=ug(2);
s(3)=xg(2);
end

if orient == 5 %climb/descend turn
s(1)=ug(1);
s(2)=ug(2);
s(3)=ug(3);
s(4)=ug(4);
s(5)=xg(2);
s(6)=xg(4);
s(7)=xg(5);
end

while yesno == 1
   options = [0 1.0E-9 1.0E-9 0 0 0 0 0 0 0 0 0 0 200000];
%    ot = input('Required# of iterations (def. = 1000): ');

%    if(isempty(ot)), options(14) = 100000;else,options(14) = ot;end
   
%    [s,~,x,u,fcost,lcost] = fminsa(cost,s,options,[],xg,ug,const,model);
   [s,~,x,u,fcost,lcost] = fminsa(cost,s,options,[],xg,ug,const,model,acft);
   [amach,qbar]=adc(x(1),x(12));
   [xd,ay,az] = model(x,u,acft);
   fprintf('\n')

   if ndof > 3
      fprintf('Throttle (percent):		%g\n', u(1))
      fprintf('Elevator (deg):			%g\n', rtod*u(2))
      fprintf('Ailerons (deg):			%g\n', rtod*u(3))
      fprintf('Rudder (deg):			%g\n', rtod*u(4))
      fprintf('Angle of Attack (deg):		%g\n', rtod*x(2))
      fprintf('Sideslip Angle (deg):		%g\n', rtod*x(3))
      fprintf('Pitch Angle (deg):		%g\n', rtod*x(5))
      fprintf('Bank Angle (deg):		%g\n', rtod*x(4))
      fprintf('Normal Acceleration:	%g\n', az)
      fprintf('Lateral Accereration:	%g\n', ay)
      fprintf('Dynamic Pressure :		%g\n', qbar)
      fprintf('Mach Number:			%g\n', amach)
      
      fprintf('Altitude:     %g\n',x(12))
   else
      fprintf('Throttle (percent):		%g\n', u(1))
      fprintf('Elevator (deg):			%g\n', u(2))
      fprintf('Alpha (deg):			    %g\n', x(2)*rtod)
      fprintf('Pitch Angle (deg):		%g\n', x(5)*rtod)
      fprintf('Normal Acceleration (g):	%g\n', az/32.2)
      fprintf('Dynamic Pressure (psf):	%g\n', qbar)
      fprintf('Mach Number:			    %g\n', amach)
   end

  fprintf('\n')
  fprintf('Initial Cost Function:		%g\n', fcost)
  fprintf('Final Cost Function:		%g\n', lcost)

%   yesno = menu('More Iterations?','Yes','No');
  yesno = 0;
   
end
Xequil=x
Uequil=u;

