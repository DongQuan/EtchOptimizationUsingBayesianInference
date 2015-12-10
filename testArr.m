function [F] = testArr(current,i)
global expParameters
k = 10e+11*(current(1)*exp(-current(2)/expParameters(i)))+lognrnd(0,current(3));
F = k;