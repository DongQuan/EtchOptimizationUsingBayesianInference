function [F]  = Posterior(current,i)
global posteriorRecord
posteriorTime = tic;
[L] = Likelihood(current);
posterior = L + Prior(current,i);
F = posterior;
posteriorTimeElapsed = toc(posteriorTime);
posteriorRecord = [posteriorRecord posterior];
assignin('base', 'posteriorTime', posteriorTimeElapsed);
end
