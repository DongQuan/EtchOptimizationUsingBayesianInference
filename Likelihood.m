function [F PlasmaParams] = Likelihood(current,expected)
global data;
global NoTrainingCases
L = 0;
for i=1:NoTrainingCases
    [sim PlasmaParams] = SimEtch(current,i,expected);
    L = log(normpdf(data(i),sim,current(15))) + L;
end
F = L;
