function [F, chainSD] = ProposalFunction(theta,current,subBlock)
global center
global proposalSD
global proposedParameterRecord
global kNorm
global proposalCenter
global proposalLB
global proposalUB
parameter = current;

parameter(subBlock) = unifrnd(proposalLB(subBlock),proposalUB(subBlock));
chainSD = 1;
%parameter(subBlock) = lognrnd(proposalCenter(subBlock),proposalSD(subBlock));
%if subBlock == 1
%    for index = 1:1
        %parameter(subBlock) = lognrnd(proposalCenter(subBlock),proposalSD(subBlock));
%         parameter(index) = parameter(index)/kNorm;
%     end
% %elseif subBlock == 2
%     for actIndex = 2:4
%         parameter(actIndex) = lognrnd(proposalCenter(actIndex),proposalSD(actIndex));
%     end
% %else
%      parameter(5) = lognrnd(proposalCenter(5),proposalSD(5));
%end
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
