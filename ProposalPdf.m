function F = ProposalPdf(current,new,subBlock,chainSD)
proposalPDFTime = tic;
%use current instead of center to change back to random walk
global proposalLambda;
global proposalSD
global proposalRecord
global proposalCenter
global proposalLB
global proposalUB
q = 1;
q = unifpdf(new(subBlock),proposalLB(subBlock),proposalUB(subBlock));
%q = lognpdf(new(subBlock),proposalCenter(subBlock),proposalSD(subBlock))*q;
% if subBlock == 1
%     for index = 1:1
%         q = lognpdf(new(index),proposalCenter(index),proposalSD(subBlock))*q;
%     end
% elseif subBlock == 2
%     for actIndex = 2:2
%         q = lognpdf(new(actIndex),proposalCenter(actIndex),proposalSD(subBlock))*q;
%     end
% else
%      q = lognpdf(new(3),proposalCenter(3),proposalSD(subBlock))*q;
% end
F = log(q);
proposalPDFTimeElapsed = toc(proposalPDFTime);
assignin('base', 'proposalPDFTimeElapsed', proposalPDFTimeElapsed);
end