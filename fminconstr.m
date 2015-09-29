function [c,ceq] = fminconstr(x,Act,B,i)
c=[];
%c = [x(9)*m/q_e*(k7_nond*n0_nond*x(11))^2-max_Te;min_Te-x(9)*m/q_e*(k7_nond*n0_nond*x(11))^2; 1 - x(13)]; %nonlinear inequality
ceq = GlobalPlasmaSystem(x,Act,B,i); % the fsolve objective is fmincon constraints