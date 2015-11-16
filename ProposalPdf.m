function F = ProposalPdf(current,new,subBlock)
proposalPDFTime = tic;
%use current instead of center to change back to random walk
global proposalLambda;
global proposalSD
global proposalRecord
global center
q = 1;
if subBlock == 1
    for index = 1:7
        q = lognpdf(new(index),center(index),proposalSD(subBlock))*q;
    end
elseif subBlock == 2
    for actIndex = 8:14
        q = lognpdf(new(actIndex),center(actIndex),proposalSD(subBlock))*q;
    end
else
     q = lognpdf(new(15),center(15),proposalSD(subBlock))*q;
end
F = q;
proposalPDFTimeElapsed = toc(proposalPDFTime);
assignin('base', 'proposalPDFTimeElapsed', proposalPDFTimeElapsed);
end