function F = ProposalFunction(current,subBlock)
global center
global proposalSD
global proposedParameterRecord
global kNorm
global proposalCenter
parameter = current;


if subBlock == 1
    for index = 1:1
        parameter(index) = lognrnd(proposalCenter(index),proposalSD(subBlock));
        parameter(index) = parameter(index)/kNorm;
    end
elseif subBlock == 2
    for actIndex = 2:2
        parameter(actIndex) = lognrnd(proposalCenter(actIndex),proposalSD(subBlock));
    end
else
     parameter(3) = lognrnd(proposalCenter(3),proposalSD(subBlock));
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
