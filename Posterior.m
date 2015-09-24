function [F ,PlasmaParams]  = Posterior(current,i,expected)
[L , PlasmaParams] = Likelihood(current,expected);
posterior = L + Prior(current,i);
F = posterior;
end
