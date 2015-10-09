global expParameters
%Generate Synthetic Data
k7_nond = 5e-8;
k1 = 3e-16;
k2 = 2.1e-18;
k3 = 1.5e-16;
k4 = 3e-15;
k5 = 1e-16;
k6 = 2e-17;
k8 = k2;
noise = 10;
expParameters = FactorialDesign();
current = [k1 k2 k3  k4 k5 k6 k8 2 6 3 2 4 1 2 0];
xmulti=zeros(25,10)
predER = zeros(1,10);
for trial=1:10
    for expNo=1:1%length(expParameters)
        xmulti(:,trial) = GlobalSolver(current,expNo);
    %PredER(trial) = CalcEtchRate(xmulti(:,trial),expNo);
    end
end
SyntheticDataWithNoise = SyntheticData + noise*randn(length(SyntheticData));


