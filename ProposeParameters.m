function F = ProposeParameters(i)
global mean
global sd
parameter = lognrnd(mean(i),sd(i));
kNorm = 10e+18;
%scale prexponential factor values back down
if i<=7
    parameter = parameter/kNorm;
end
%parameter = abs(parameter);
F = parameter;
end
