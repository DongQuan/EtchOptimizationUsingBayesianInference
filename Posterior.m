function [F]  = Posterior(current,i)
[L] = Likelihood(current);
posterior = L + Prior(current,i);
F = posterior;
end
