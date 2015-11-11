%Test function for Metropolis Hastings Method
clear 
close all
pause('on')
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
global noExpParameters
global noUnknowns
global center
global sd
global proposalLambda
global priorRecord
global proposalRecord
global posteriorRecord
global alphaRecord
global proposedParameterRecord
global likelihoodRecord
global etchRecord

%define constants

massIon = 35.45*(1.660568e-27); %kg
kb = 1.381e-23; %J/K
q = 1.6022e-19; %C
R = 0.15; %m
L = 0.14; %m
T = 600;
V = pi*R^2*L;
A = 2*pi*R*L + 2*pi*R^2;
K = 5e-14; %(m^3/s)
plasmaUnknowns = 25;
sigma = 16.8e-24;
noUnknowns = 15; %7 k's (A and B coeff) plus standard error
proposalLambda = [.1 .1 .1 .1 .1 .1 .1 .01 .01 .01 .01 .01 .01 .01 .01 1];
priorRecord = [];
proposalRecord = [];
posteriorRecord = [];
alphaRecord = [];
proposedParameterRecord = [];
likelihoodRecord = [];
etchRecord = [];

k7_nond = 5e-8;
k1 = log(500);
k2 = log(5);
k3 = log(500);
k4 = log(5000);
k5 = log(5);
k6 = log(50);
k8 = k2;
center= [k1 k2 k3  k4 k5 k6 k8 log(5) log(5) log(5) log(5) log(5) log(5) log(5) log(10)];
sd = [1 1 1 1 1 1 1 .7 .7 .7 .7 .7 .7 .7 1];
noExpParameters = 7;
subBlocks = 3;
%Calculate k9
diffusionLength = sqrt(1/(2.405/R)^2+(pi/L)^2);
Df = diffusionLength/3*sqrt(8*kb*T/(pi*massIon)); %need to calculate effective diffusion coefficient
recombinationProbability = (0.2+10e-3+0.05)/3;
k9 = recombinationProbability*Df/(diffusionLength^2);

%Split data into training and test sets
allExpParameters = xlsread('SyntheticData.xlsx','ExpParameters');
trainingDataIndex = randperm(size(allEtchRates,1)+1);
trainingDataIndex = trainingDataIndex(1:size(allEtchRates,1)/6)-1;
noTrainingCases = length(trainingDataIndex);
noTestCases = size(allEtchRates,1) - noTrainingCases;
trainingData = zeros(noTrainingCases,1);
trainingExp = zeros(noTrainingCases,noExpParameters);
testData = [];
testExp = [];
for i=1:noTrainingCases
    trainingData(i) = allEtchRates(trainingDataIndex(i));
    trainingExp(i,:) = allExpParameters(trainingDataIndex(i),:);
end
Lia = ismember(allEtchRates,trainingData);
for i=1:length(Lia)
    if ~(Lia(i))
        testData = [testData allEtchRates(i)];
        testExp = [testExp; allExpParameters(i,:)];
    end
end
%assign training data and exp parameters to global variables
%data = trainingData;%cat(1,trainingData,testData);
% expParameters = trainingExp; %cat(1,trainingExp,testExp);

%Make initial guess
current = zeros(noUnknowns,1);
for i = 1:subBlocks
      current = ProposeParameters(current,i);
end
%current(1) = 7e-12;
%current(2) = 7;
%current = [2.34E-16	1.13E-19	1.73E-17	2.70E-16	4.12E-19	2.33E-18	1.21E-19	5.502247949	5.297391123	3.862552087	4.102946538	54.31011931	8.675632081	3.320856533	7.35781683];
%set experiments to training data


data = realValues;
index = 1;
nn      = 100;       % Number of samples for examine the AC
N       = 100;     % Number of samples (iterations)
burnin  = 1;      % Number of runs until the chain approaches stationarity
lag     = 1;        % Thinning or lag period: storing only every lag-th point
theta   = zeros(N*subBlocks,noUnknowns); 
acc = zeros(subBlocks,1);

%initialize plasma arrays
AlphaSet = zeros(N,subBlocks);
% EtchRateSetMean = zeros(1, length(trainingData));
% TestEtchRateSetMean = zeros(1, length(testData));
% EtchRateSetMode = zeros (1,length(trainingData));
% for i = 1:burnin    % First make the burn-in stage
%     [t] = MetropolisHastings(current);
% end
count = 0;
[PosteriorCurrent] = Posterior(current,1);
%Peform MH iterations
totalTime = tic;

for cycle = 1:N  % Cycle to the number of samples
    %for j = 1:lag 
    MHtime = tic;
    for j=1:2 % Cycle to make the thinning
        SCtime = tic;
        [alpha,t, a,prob, PosteriorCatch] = MetropolisHastings(current,PosteriorCurrent,j);
        SCelapsed = toc(SCtime);
        theta(index,:) = t;        % Samples accepted
        AlphaSet(cycle,j) = alpha;
        index = index + 1;
        current = t;
        PosteriorCurrent = PosteriorCatch;
        acc(j) = acc(j) + a;  % Accepted ?
    end
    count = count+1;
    MHelapsed = toc(MHtime);
end
totalTimElapsed = toc(totalTime);
accrate = acc/N;     % Acceptance rate
% hf2 = mode(theta);
% hfmean = zeros(1,NoUnknowns);
% for i = 1:NoUnknowns
%     hfmean(i) = sum(theta(:,i))/N;

simTrainingData = zeros(length(trainingData),1);
for i=1:length(data)
    [simTrainingData(i),simPlasmaParams(:,i)] = GlobalSolver(mean(theta),i);
end
figure(1)
x = linspace(1,length(trainingData),length(trainingData));
scatter(x,simTrainingData,'r')
hold on
scatter(x,realValues,'g')
expParameters = testExp;
for i=1:length(testExp)
    [realTestValues(i),realtestPlasmaParams(:,i)]= GlobalSolver(real,i);
    [simTestData(i),simTestPlasmaParams(:,i)] = GlobalSolver(mean(theta),i);
end

figure(2)
x = linspace(1,length(testData),length(testData));
scatter(x,simTestData,'r')
hold on
scatter(x,realTestValues,'g')

