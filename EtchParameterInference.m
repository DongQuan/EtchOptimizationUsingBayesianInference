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
mean = [k1 k2 k3  k4 k5 k6 k8 2 6 3 2 4 1 2 0];
sd(15) = 0;
allExpParameters = FactorialDesign();
noExpParameters = 7;

%Calculate k9
diffusionLength = sqrt(1/(2.405/R)^2+(pi/L)^2)
Df = diffusionLength/3*sqrt(8*kb*T/(pi*massIon)); %need to calculate effective diffusion coefficient
recombinationProbability = (0.2+10e-3+0.05)/3;
k9 = recombinationProbability*Df/(diffusionLength^2);

%Split data into training and test sets
AllData = LoadSynData();
TrainingDataIndex = randperm(size(allExpParameters,1)+1);
TrainingDataIndex = TrainingDataIndex(1:size(allExpParameters,1)/2)-1;
NoTrainingCases = length(TrainingDataIndex);
NoTestCases = size(allExpParameters,1) - NoTrainingCases;
TrainingData = zeros(NoTrainingCases,1);
TrainingExp = zeros(NoTrainingCases,noExpParameters);
TestData = zeros(NoTestCases,1);
TestExp = zeros(NoTestCases,noExpParameters);
for i=1:NoTrainingCases
    TrainingData(i) = AllData(TrainingDataIndex(i));
    TrainingExp(i,:) = allExpParameters(TrainingDataIndex(i),:);
end
for i=1:NoTestCases
    if i~=(TrainingDataIndex(i))
    TestData(i) = AllData(i);
    TestExp(i,:) = allExpParameters(i,:);
    end
end
%Merge training and test data to one matrix
data = cat(1,TrainingData,TestData);
expParameters = cat(1,TrainingExp,TestExp);

%Make initial guess
current = zeros(noUnknowns,1);

for i = 1:noUnknowns
     current(i) = ProposeParameters(i);
end


[PosteriorCurrent] = Posterior(current,1);

index = 1;
nn      = 100;       % Number of samples for examine the AC
N       = 1;     % Number of samples (iterations)
burnin  = 1;      % Number of runs until the chain approaches stationarity
lag     = 1;        % Thinning or lag period: storing only every lag-th point
theta   = zeros(N*noUnknowns,noUnknowns); 
acc = 0;

%initialize plasma arrays
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

% %Peform MH iterations
% for cycle = 1:N   % Cycle to the number of samples
%     %for j = 1:lag 
%     for j=1:1%NoUnknowns % Cycle to make the thinning
%         [alpha,t, a,prob, PosteriorCatch] = MetropolisHastings(current,PosteriorCurrent,j,expected);
%         theta(index,:) = t;        % Samples accepted
%         AlphaSet(i,j) = alpha;
%         index = index + 1;
%         current = t;
%         PosteriorCurrent = PosteriorCatch;
%     end
%     PlasmaParamsSet(:,i) = expected;
%     acc      = acc + a;  % Accepted ?
%     count = count+1
% end
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
