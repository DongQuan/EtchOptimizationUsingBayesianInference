function [F] = testArr(current,i)
global expParameters
k = (current(1)*exp(-current(2)/expParameters(i,2)))+current(3)^2 + current(4)*sin(current(6)*pi*expParameters(i,2))+ sqrt(current(5)) + normrnd(0,current(7));
F = k;

