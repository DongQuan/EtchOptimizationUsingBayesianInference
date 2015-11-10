function F = ProposalPdf(current,new,i)
proposalPDFTime = tic;
prior = 0;
%use current instead of mean to change back to random walk
global proposalLambda;
global center
global proposalRecord
q = log(exppdf(proposalLambda(i),current(i))+10e-20);
proposalRecord = [proposalRecord q];
F = q;
proposalPDFTimeElapsed = toc(proposalPDFTime);
assignin('base', 'proposalPDFTimeElapsed', proposalPDFTimeElapsed);
end