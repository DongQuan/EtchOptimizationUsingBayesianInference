function [alpha, t, a, prob,PosteriorCatch] = MetropolisHastings(current,PosteriorCurrent,i)
new = current;
new(i) = ProposeParameters(i);
[PosteriorNew] = Posterior(new,i)
% sampling from the proposal PDF with media the current state
[alpha] = exp(PosteriorNew + ProposalPdf(current,new,i)-(PosteriorCurrent+ProposalPdf(new,current,i))) % Ratio of the density at the candidate (theta_ast) and current (current) points
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
