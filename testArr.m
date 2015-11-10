function [F] = testArr(current)
k = 10e+12*current(1)*exp(-current(2));
F = k;