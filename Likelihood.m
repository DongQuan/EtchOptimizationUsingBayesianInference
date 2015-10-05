function [F] = Likelihood(current)
global data;
L = 0;
for i=1:length(data)
    [sim] = GlobalSolver(current,i);
    L = log(normpdf(data(i),sim,current(15))) + L;
end
F = L;
