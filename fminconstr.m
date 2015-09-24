function [c,ceq] = fminconstr(x,A,B,i)
global ExpParameters
global m
global kb
global q_e
global T
P0_nond = ExpParameters(i,1);
T0_nond = T;
n0_nond = P0_nond/(kb*T0_nond)/(100^3);%convert to cm^3
k7_nond = 5e-8;
max_Te = 5;
min_Te = 0;
%c=[];
c = [x(9)*m/q_e*(k7_nond*n0_nond*x(11))^2-max_Te;min_Te-x(9)*m/q_e*(k7_nond*n0_nond*x(11))^2; 1 - x(13)]; %nonlinear inequality
ceq = GlobalPlasmaSystem(x,A,B,i); % the fsolve objective is fmincon constraints