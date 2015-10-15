function [c,ceq] = fminconstr(x,Act,B)
c=[];
ceq = GlobalPlasmaSystemWithDimensionsAndK(x,Act,B); % the fsolve objective is fmincon constraints