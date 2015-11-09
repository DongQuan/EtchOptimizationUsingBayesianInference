function F = ProposalPdf(current,new,i)
proposalPDFTime = tic;
prior = 0;
%use current instead of mean to change back to random walk
global proposalSigma;
global center
global proposalRecord
hold = 77
[current(i)]
[new(i)]

kNorm = 10e+18;
% if i<=7
%     mu = mu/kNorm;
%     sigma = sigma/kNorm;
% end

q = log(lognpdf(new(i),current(i),proposalSigma(i))+10e-20)
hold = 77
proposalRecord = [proposalRecord q];
F = q;
proposalPDFTimeElapsed = toc(proposalPDFTime);
assignin('base', 'proposalPDFTimeElapsed', proposalPDFTimeElapsed);
end