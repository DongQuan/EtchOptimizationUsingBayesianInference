function F = ProposalPdf(a,b,i)
prior = 0;
global sd;
q = log(normpdf(a,b,sd(i))) + prior;
F = q;
end