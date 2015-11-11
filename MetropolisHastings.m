function [alpha, t, a, prob,PosteriorCatch] = MetropolisHastings(current,PosteriorCurrent,subBlock)
new = ProposeParameters(current,subBlock);
[PosteriorNew] = Posterior(new,subBlock)
% current
% new
% prior1 = Prior(current,subBlock)
% prior2 = Prior(new,subBlock)
% p1 = ProposalPdf(current,new,subBlock)
% l1 = Likelihood(current)
% l2 = Likelihood(new)
% p2 = ProposalPdf(new,current,subBlock)
% sampling from the proposal PDF with media the current state
[alpha] = exp(PosteriorNew + ProposalPdf(current,new,subBlock)-(PosteriorCurrent+ProposalPdf(new,current,subBlock))) % Ratio of the density at the candidate (theta_ast) and current (current) points
%alphaRecord = [alphaRecord alpha];
if rand <= min(alpha,1)
   t    = new;        % Accept the candidate
   prob = min(alpha,1);     % Accept with probability min(alpha,1)
   a    = 1;                % Note the acceptance
   PosteriorCatch=PosteriorNew;
else
   a = 0;
   t    = current;            % Reject the candidate and use the same state
   prob = 1-min(alpha,1);   % The same state with probability 1-min(alpha,1)
   PosteriorCatch = PosteriorCurrent;
end
end
