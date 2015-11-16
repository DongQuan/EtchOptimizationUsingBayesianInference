function F = ProposalFunction(current,subBlock)
global center
global proposalSD
global proposedParameterRecord
parameter = current;
kNorm = 10e+18;

if subBlock == 1
    for index = 1:7
        parameter(index) = lognrnd(center(index),proposalSD(subBlock));
        parameter(index) = parameter(index)/kNorm;
    end
elseif subBlock == 2
    for actIndex = 8:14
        parameter(actIndex) = lognrnd(center(actIndex),proposalSD(subBlock));
    end
else
     parameter(15) = lognrnd(center(15),proposalSD(subBlock));
end
% if subBlock == 1
%     for (index = 1:7)
%      parameter(index)= unifrnd(0,10000)/kNorm;
%     end
% elseif subBlock == 2
%      for (index = 8:14)
%         parameter(index)= unifrnd(0,30);
%      end
% else
%     parameter(15)= unifrnd(0,50);
% end

proposedParameterRecord = [proposedParameterRecord parameter];
F = parameter;
end
