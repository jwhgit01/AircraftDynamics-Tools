% Controllability and Observability
format shortg
close all
clear
clc
addpath(genpath('src/'))

% LTV system
n = 2;
m = 1;
p = 2;
A = @(t) [-2, 1+sin(t); 0, -4];
B = @(t) [0;cos(t)];
C = @(t) [1, 0; 0, 1+cos(t)];
D = [0;0];
T = 30;

% Controllability Gramian
[tC,WC] = ControllabilityGramian(A,B,[0 T]);
NC = length(tC);
RelambdaC = zeros(NC,n);
for ii = 1:NC
    RelambdaC(ii,:) = real(eig(squeeze(WC(ii,:,:)))).';
end
figure
plot(tC,RelambdaC)
ylabel('$\mathrm{Re}(\lambda(W_C))$','Interpreter','latex')
xlabel('Time')
title('Controllability Gramian')

% Observability Gramian
[tO,WO] = ObservabilityGramian(A,C,[0 T]);
NO = length(tO);
RelambdaO = zeros(NO,n);
for ii = 1:NO
    RelambdaO(ii,:) = real(eig(squeeze(WO(ii,:,:)))).';
end
figure
plot(tO,RelambdaO)
ylabel('$\mathrm{Re}(\lambda(W_O))$','Interpreter','latex')
xlabel('Time')
title('Observability Gramian')