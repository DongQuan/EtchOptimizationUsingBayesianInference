function F = ProposeParameters(new,subBlock)
global proposalCenter
global proposalSD
global proposedParameterRecord
global kNorm
parameter = new;

if subBlock == 1
    for index = 1:1
        parameter(index) = lognrnd(proposalCenter(index),proposalSD(index));
        parameter(index) = parameter(index)/kNorm;
    end
elseif subBlock == 2
    for actIndex = 2:2
        parameter(actIndex) = lognrnd(proposalCenter(actIndex),proposalSD(actIndex));
    end
else
     parameter(3) = lognrnd(proposalCenter(3),proposalSD(3));
end

proposedParameterRecord = [proposedParameterRecord parameter];
F = parameter;
end
