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
global noExpParameters
global noUnknowns
global mean
global sd


%define constants

massIon = 35.45*(1.660568e-27); %kg
kb = 1.381e-23; %J/K
q = 1.6022e-19; %C
R = 0.15; %m
L = 0.14; %m
T = 303;
V = pi*R^2*L;
A = 2*pi*R*L + 2*pi*R^2;
K = 5e-14; %(m^3/s)
plasmaUnknowns = 25;
sigma = 16.8e-24;
noUnknowns = 15; %7 k's (A and B coeff) plus standard error

sd = zeros(noUnknowns,1);
k7_nond = 5e-8;
k1 = 3e-16;
k2 = 2.1e-18;
k3 = 1.5e-16;
k4 = 3e-15;
k5 = 1e-16;
k6 = 2e-17;
k8 = k2;
mean = [k1 k2 k3  k4 k5 k6 k8 2 6 3 2 4 1 2 20];
sd = [1 1 1 1 1 1 1 5 5 5 5 5 5 5 20];
noExpParameters = 8;

%Calculate k9
diffusionLength = sqrt(1/(2.405/R)^2+(pi/L)^2)
Df = diffusionLength/3*sqrt(8*kb*T/(pi*massIon)); %need to calculate effective diffusion coefficient
recombinationProbability = (0.2+10e-3+0.05)/3;
k9 = recombinationProbability*Df/(diffusionLength^2);

%Split data into training and test sets
allEtchRates = xlsread('SyntheticData.xlsx','EtchRates');
allExpParameters = xlsread('SyntheticData.xlsx','ExpParameters');
trainingDataIndex = randperm(size(allEtchRates,1)+1);
trainingDataIndex = trainingDataIndex(1:size(allEtchRates,1)/2)-1;
noTrainingCases = length(trainingDataIndex);
noTestCases = size(allEtchRates,1) - noTrainingCases;
trainingData = zeros(noTrainingCases,1);
trainingExp = zeros(noTrainingCases,noExpParameters);
testData = zeros(noTestCases,1);
testExp = zeros(noTestCases,noExpParameters);
for i=1:noTrainingCases
    trainingData(i) = allEtchRates(trainingDataIndex(i));
    trainingExp(i,:) = allExpParameters(trainingDataIndex(i),:);
end
for i=1:noTestCases
    if i~=(trainingDataIndex(i))
    testData(i) = allEtchRates(i);
    testExp(i,:) = allExpParameters(i,:);
    end
end
%assign training data and exp parameters to global variables
data = trainingData;%cat(1,trainingData,testData);
expParameters = trainingExp; %cat(1,trainingExp,testExp);

%Make initial guess
current = zeros(noUnknowns,1);

for i = 1:noUnknowns
     current(i) = ProposeParameters(i);
end


[PosteriorCurrent] = Posterior(current,1);

index = 1;
nn      = 100;       % Number of samples for examine the AC
N       = 300;     % Number of samples (iterations)
burnin  = 1;      % Number of runs until the chain approaches stationarity
lag     = 1;        % Thinning or lag period: storing only every lag-th point
theta   = zeros(N*noUnknowns,noUnknowns); 
acc = 0;

%initialize plasma arrays
AlphaSet = zeros(N,noUnknowns);
EtchRateSetMean = zeros(1, length(trainingData));
TestEtchRateSetMean = zeros(1, length(testData));
EtchRateSetMode = zeros (1,length(trainingData));
% for i = 1:burnin    % First make the burn-in stage
%     [t] = MetropolisHastings(current);
% end
count = 0;

%Peform MH iterations
for cycle = 1:N   % Cycle to the number of samples
    %for j = 1:lag 
    for j=1:noUnknowns % Cycle to make the thinning
        [alpha,t, a,prob, PosteriorCatch] = MetropolisHastings(current,PosteriorCurrent,j);
        theta(index,:) = t;        % Samples accepted
        AlphaSet(i,j) = alpha;
        index = index + 1;
        current = t;
        PosteriorCurrent = PosteriorCatch;
    end
    acc      = acc + a;  % Accepted ?
    count = count+1
end
% accrate = acc/N;     % Acceptance rate
% hf2 = mode(theta);
% hfmean = zeros(1,NoUnknowns);
% for i = 1:NoUnknowns
%     hfmean(i) = sum(theta(:,i))/N;
% end
% for j = 1:length(TrainingData)
%          EtchRateSetMean(j) = SimEtch(hfmean,j,expected);
% end
% ind = 1;
% for j = length(TrainingData)+1:length(TrainingData)+length(TestData)
%          TestEtchRateSetMean(ind) = SimEtch(hfmean,j,expected);
%          ind = ind +1;
% end
% figure(1)
% x = linspace(1,length(TrainingData),length(TrainingData));
% scatter(x,EtchRateSetMean)
% hold on
% scatter(x,TrainingData,'g')
% print('-f1','EtchSimulation','-dpng')
% fid = fopen('theta.txt', 'wt'); % Open for writing
% for i=1:size(theta,1)
%    fprintf(fid, '%d ', theta(i,:));
%    fprintf(fid, '\n');
% end
% fclose(fid);
% 
% %Test set
% figure(2)
% x = linspace(1,length(TestData),length(TestData));
% scatter(x,TestEtchRateSetMean)
% hold on
% scatter(x,TestData,'g')
% print('-f1','EtchPrediction','-dpng')
% fid = fopen('theta.txt', 'wt'); % Open for writing
% for i=1:size(theta,1)
%    fprintf(fid, '%d ', theta(i,:));
%    fprintf(fid, '\n');
% end
%fclose(fid);
