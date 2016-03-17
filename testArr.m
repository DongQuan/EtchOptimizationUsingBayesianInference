function [F] = testArr(current,i)
global expParameters
k = (current(1)*exp(-current(2)/expParameters(i)))+normrnd(0,current(3));
F = k;