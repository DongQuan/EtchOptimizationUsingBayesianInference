function F = ProposalPdf(current,new,subBlock)
proposalPDFTime = tic;
%use current instead of mean to change back to random walk
global proposalLambda;
global center
global proposalRecord
q = 0;
sigma1 = eye(7)*4;
sigma2 = eye(7)*2;
sigma3 = eye;
if subBlock==1
    q = mvnpdf(current(1:7),new(1:7),sigma1);
elseif subBlock==2
    q = mvnpdf(current(8:14),new(8:14),sigma2);
else 
    q = mvnpdf(current(15),new(15),sigma3);
end
proposalRecord = [proposalRecord q];
F = q;
proposalPDFTimeElapsed = toc(proposalPDFTime);
assignin('base', 'proposalPDFTimeElapsed', proposalPDFTimeElapsed);
end