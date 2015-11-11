function F = ProposeParameters(new,subBlock)
global center
global sd
global proposedParameterRecord
parameter = new;
kNorm = 10e+18;
if subBlock == 1
    for index = 1:7
        parameter(index) = lognrnd(center(index),sd(index));
        parameter(index) = parameter(index)/kNorm;
    end
elseif subBlock == 2
    for actIndex = 8:14
        parameter(actIndex) = lognrnd(center(actIndex),sd(actIndex));
    end
else
     parameter(15) = lognrnd(center(15),sd(15));
end

proposedParameterRecord = [proposedParameterRecord parameter];
F = parameter;
end
