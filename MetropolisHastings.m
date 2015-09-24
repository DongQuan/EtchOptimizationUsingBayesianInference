function [alpha, t, a, prob,expected2,PosteriorCatch] = MetropolisHastings(current,PosteriorCurrent,i,expected)
new = current;
new(i) = ProposeParameters(i);   
[PosteriorNew expected_new] = Posterior(new,i,expected);
% sampling from the proposal PDF with media the current state
alpha = exp(PosteriorNew + ProposalPdf(current(i),new(i),i)-(PosteriorCurrent+ProposalPdf(new(i),current(i),i)));  % Ratio of the density at the candidate (theta_ast) and current (current) points
if rand <= min(alpha,1)
   t    = new;        % Accept the candidate
   prob = min(alpha,1);     % Accept with probability min(alpha,1)
   a    = 1;                % Note the acceptance
   PosteriorCatch=PosteriorNew;
   expected2 = expected_new
else
   a = 0;
   t    = current;            % Reject the candidate and use the same state
   prob = 1-min(alpha,1);   % The same state with probability 1-min(alpha,1)
   PosteriorCatch = PosteriorCurrent;
   expected2 = expected;
end
end
