function [c,ceq] = fminconstr(x,Act,B)
c=[];
ceq = GlobalPlasmaSystemWithDimensions(x,Act,B); % the fsolve objective is fmincon constraints