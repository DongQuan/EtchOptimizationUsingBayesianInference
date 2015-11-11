function [F]  = Posterior(current,subBlock)
global posteriorRecord
posteriorTime = tic;
[L] = Likelihood(current);
posterior = L + Prior(current,subBlock);
F = posterior;
posteriorTimeElapsed = toc(posteriorTime);
posteriorRecord = [posteriorRecord posterior];
assignin('base', 'posteriorTime', posteriorTimeElapsed);
end
