d = 5;
N = 1;
randomDraws = unidrnd(N,d,1);
global expParameters
expParameters = FactorialDesign();
current = [k1 k2 k3 k4 k5 k6 k8 0 0 0 0 0 0 0 0]
data = 1;

for j = 1:length(data)
%    for i=1:d
        trainingPred(j) = GlobalSolver(current,j);
%    end
end

k7_nond = 5e-8;
k1 = 3e-16;
k2 = 2.1e-18;
k3 = 1.5e-16;
k4 = 3e-15;
k5 = 1e-16;
k6 = 2e-17;
k8 = k2;
