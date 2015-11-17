function [alpha, t, a, prob,PosteriorCatch] = MetropolisHastings(current,PosteriorCurrent,subBlock)
%independent sampling from proposal function (current is put here to hold
%unused blocks of parameters)
new = ProposalFunction(current,subBlock);
[PosteriorNew] = Posterior(new,subBlock);
[alpha] = exp(PosteriorNew + ProposalPdf(current,new,subBlock)-(PosteriorCurrent+ProposalPdf(new,current,subBlock))); % Ratio of the density at the candidate (theta_ast) and current (current) points
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
