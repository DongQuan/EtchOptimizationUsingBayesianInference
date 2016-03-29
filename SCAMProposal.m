function [parameter,chainSD] = SCAMProposal(theta,current,index)
global proposalSD
global proposedParameterRecord
global kNorm
global proposalCenter
parameter = current;
s = 2.4;
epsilon = 10e-3;
variance = s*(var(theta(:,index)))+s*epsilon;
chainSD = sqrt(log(variance/((proposalCenter(index))^2)+1));
parameter(index) = lognrnd(proposalCenter(index),chainSD);
proposedParameterRecord = [proposedParameterRecord parameter];
end