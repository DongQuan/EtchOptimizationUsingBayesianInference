function F = ProposalPdf(current,new,i)
proposalPDFTime = tic;
prior = 0;
%use current instead of mean to change back to random walk
global sd;
mu = mean(i);
sigma = sd(i);
kNorm = 10e+18;
if i<=7
    mu = mu/kNorm;
    sigma = sigma/kNorm;
end

q = log(normpdf(new,current,sigma)) + prior;
F = q;
proposalPDFTimeElapsed = toc(proposalPDFTime);
assignin('base', 'proposalPDFTimeElapsed', proposalPDFTimeElapsed);
end