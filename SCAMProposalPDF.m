function F = SCAMProposalPDF(current,new,index,chainSD)
global center
global proposalSD
global proposedParameterRecord
global kNorm
global proposalCenter
q=1;
q = lognpdf(new(index),proposalCenter(index),chainSD)*q;
F = log(q);
end