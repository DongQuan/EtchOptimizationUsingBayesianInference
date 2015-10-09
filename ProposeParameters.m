function F = ProposeParameters(i)
global mean
global sd
parameter = normrnd(mean(i),sd(i));
parameter = abs(parameter);
F = parameter;
end
