%Test function for Metropolis Hastings Method
clear 
close all
rng('default');
% declare global variables
global massIon %ion mass
global kb %boltzmann constant
global q %electron charge
global R %radius of chamber
global L %length of chamber
global T %temperature of system (treat constant)
global V %volume
global A %area
global K %this is equal to k7 in Efremov study
global plasmaUnknowns
global k9
global sigma
global expParameters
global data




kb = 1.381*10^-23;
T = 303;
q_e = 1.6e-19;
m = 35.45*1.660468e-27; 
R = 15; 
L = 14;
NoUnknowns = 15; %7 k's (A and B coeff) plus standard error
count = 0;
sd = zeros(NoUnknowns,1);
k7_nond = 5e-8;
k1 = 3e-10/k7_nond;
k2 = 2.1e-12/k7_nond;
k3 = 1.5e-10/k7_nond;
k4 = 3e-9/k7_nond;
k5 = 1e-10/k7_nond;
k6 = 2e-11/k7_nond;
k8 = k2;
mean = [k1 k2 k3  k4 k5 k6 k8 2 6 3 2 4 1 2];
mean(15) = 0;
sd(15) = 0;
ExpParameters = FactorialDesign();
%TeSet = [3.65 3.49 3.25 3.2 3.05 3 2.95 2.9 2.85 2.83 2.8 3.25 2.75 2.5 2.4 2.3 2.25 2.75 2.6 2.4 2.23 2.78 2.65 2.55 2.45 2.3 2.26 2.24 2.23 2.21 2.2 2.2];
%ExpParameters = [1.33 303 20 700 0 0 300;1.33 303 20 700 0 1 300; 1.33 303 20 800 0 .7 250;.33 303 20 800 0 .7 300; 1.33 303 20 700 0 .5 300; 1.33 303 20 700 0 .7 300];
SynData = zeros(size(ExpParameters,1),1);
PlasmaSet = zeros(PlasmaUnknowns,length(SynData));
for ExpNo = 1:size(ExpParameters,1)
     SynData(ExpNo) = SimEtch(mean,ExpNo,expected);
end

x = linspace(1,length(SynData),length(SynData));
scatter(x,SynData)



AllData = [140.230972761902;130.410857280232;115.212188610325;112.233750440762;102.616977408830;99.5317573706136;98.2774194559154;92.6430362142662;89.3688919593096;87.8189438757718;85.9200553557556;115.212188610325;129.948827258397;102.711297162993;91.4399518046163;79.9332744171053;74.0964936062844;129.984600731042;113.783612812632;91.4396836679203;71.7117595603738;133.171616588336;119.234844964471;108.280487923322;97.1100657415797;79.9217456486127;75.2457181848829;72.8925760918235;71.7121746856209;69.3436933825009;68.1555817666271;68.1555724539584];
TrainingDataIndex = randperm(size(AllExpParameters,1)+1);
TrainingDataIndex = TrainingDataIndex(1:size(AllExpParameters,1)/2)-1;
NoTrainingCases = length(TrainingDataIndex);
NoTestCases = size(AllExpParameters,1) - NoTrainingCases;
TrainingData = zeros(NoTrainingCases,1);
TrainingExp = zeros(NoTrainingCases,NoExpParameters);
TestData = zeros(NoTestCases,1);
TestExp = zeros(NoTestCases,NoExpParameters);
for i=1:NoTrainingCases
    TrainingData(i) = AllData(TrainingDataIndex(i));
    TrainingExp(i,:) = AllExpParameters(TrainingDataIndex(i),:);
end
for i=1:NoTestCases
    if i~=(TrainingDataIndex(i))
    TestData(i) = AllData(i);
    TestExp(i,:) = AllExpParameters(i,:);
    end
end
%Merge training and test data to one matrix
data = cat(1,TrainingData,TestData);
ExpParameters = cat(1,TrainingExp,TestExp);
for i = 1:NoUnknowns
     current(i) = ProposeParameters(i);
end

%Test new system of equations
for j = 1:length(TrainingData)
         EtchRateSetMean(j) = SimEtch(mean,j,expected);
end
[PosteriorCurrent] = Posterior(current,1,expected);
% Parameters
index = 1;
nn      = 100;       % Number of samples for examine the AC
N       = 1;     % Number of samples (iterations)
burnin  = 1;      % Number of runs until the chain approaches stationarity
lag     = 1;        % Thinning or lag period: storing only every lag-th point
theta   = zeros(N*NoUnknowns,NoUnknowns); 
acc = 0;
AlphaSet = zeros(N,NoUnknowns);
EtchRateSetMean = zeros(1, length(TrainingData));
TestEtchRateSetMean = zeros(1, length(TestData));
EtchRateSetMode = zeros (1,length(TrainingData));
DimensionlessSet = zeros(17,N);
PlasmaParamsSet=zeros(PlasmaUnknowns,N);
% for i = 1:burnin    % First make the burn-in stage
%     [t] = MetropolisHastings(current);
% end
count = 0;

for cycle = 1:N   % Cycle to the number of samples
    %for j = 1:lag 
    for j=1:1%NoUnknowns % Cycle to make the thinning
        [alpha,t, a,prob, expected,PosteriorCatch] = MetropolisHastings(current,PosteriorCurrent,j,expected);
        theta(index,:) = t;        % Samples accepted
        AlphaSet(i,j) = alpha;
        index = index + 1;
        current = t;
        PosteriorCurrent = PosteriorCatch;
    end
    PlasmaParamsSet(:,i) = expected;
    acc      = acc + a;  % Accepted ?
    count = count+1
end
accrate = acc/N;     % Acceptance rate
hf2 = mode(theta);
hfmean = zeros(1,NoUnknowns);
for i = 1:NoUnknowns
    hfmean(i) = sum(theta(:,i))/N;
end
for j = 1:length(TrainingData)
         EtchRateSetMean(j) = SimEtch(hfmean,j,expected);
end
ind = 1;
for j = length(TrainingData)+1:length(TrainingData)+length(TestData)
         TestEtchRateSetMean(ind) = SimEtch(hfmean,j,expected);
         ind = ind +1;
end
figure(1)
x = linspace(1,length(TrainingData),length(TrainingData));
scatter(x,EtchRateSetMean)
hold on
scatter(x,TrainingData,'g')
print('-f1','EtchSimulation','-dpng')
fid = fopen('theta.txt', 'wt'); % Open for writing
for i=1:size(theta,1)
   fprintf(fid, '%d ', theta(i,:));
   fprintf(fid, '\n');
end
fclose(fid);

%Test set
figure(2)
x = linspace(1,length(TestData),length(TestData));
scatter(x,TestEtchRateSetMean)
hold on
scatter(x,TestData,'g')
print('-f1','EtchPrediction','-dpng')
fid = fopen('theta.txt', 'wt'); % Open for writing
for i=1:size(theta,1)
   fprintf(fid, '%d ', theta(i,:));
   fprintf(fid, '\n');
end
fclose(fid);
