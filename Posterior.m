function [F]  = Posterior(current,i)
posteriorTime = tic;
[L] = Likelihood(current);
posterior = L + Prior(current,i);
F = posterior;
posteriorTimeElapsed = toc(posteriorTime);
assignin('base', 'posteriorTime', posteriorTimeElapsed);
end
