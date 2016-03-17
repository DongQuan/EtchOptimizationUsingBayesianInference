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
global proposalSD

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
noUnknowns = 3; %7 k's (A and B coeff) plus standard error
noExpParameters = 7;
global proposalCenter;
global kNorm

kNorm = 1;%10e+10;
subBlocks = 3;
%Calculate k9
diffusionLength = sqrt(1/(2.405/R)^2+(pi/L)^2);
Df = diffusionLength/3*sqrt(8*kb*T/(pi*massIon)); %need to calculate effective diffusion coefficient
recombinationProbability = (0.2+10e-3+0.05)/3;
k9 = recombinationProbability*Df/(diffusionLength^2);

%initialize global variables for checking Metropolis Hastings
priorRecord = [];
proposalRecord = [];
posteriorRecord = [];
alphaRecord = [];
proposedParameterRecord = [];
likelihoodRecord = [];
etchRecord = [];

%Test Arrhenius parameters
testr1 = 13;
testr2 = 5;
error = .2;
real = [testr1 testr2 error];
center = [10 5 .2];
sd = [2 2 .2]
%sd = ones(15,1)*.25;
for i = 1:length(center)
    m = center(i);
    v = sd(i)^2;
    proposalSD(i) = sqrt(log(v/(m^2)+1));
    proposalCenter(i) = log((m^2)/sqrt(v+m^2));
end
%X = lognrnd(proposalCenter(2),proposalSD(2),1,30)
% %Define real values for unknowns to generate syntehtic data
% r1 = 3e-16;
% r2 = 2.1e-18;
% r3 = 1.5e-16;
% r4 = 3e-15;
% r5 = 1e-16;
% r6 = 2e-17;
% r8 = r2;
% expError = 10;
% real = [r1 r2 r3  r4 r5 r6 r8 2 5 6 10 18 5 8 5];
% 
% %initial guesses for Metropolis Hastings
% k1 = log(500);
% k2 = log(5);
% k3 = log(500);
% k4 = log(5000);
% k5 = log(5);
% k6 = log(50);
% k8 = k2;
% 
% %Priors
% center= [k1 k2 k3  k4 k5 k6 k8 log(5) log(5) log(5) log(5) log(5) log(5) log(5) log(10)];
% sd = ones(15);
% sd(15) = 5;
% 
% load synthetic
% load allExpParameters
% %sd = [1 1 1 1 1 1 1 .7 .7 .7 .7 .7 .7 .7 1
% for i=1:length(synthetic)
% synthetic(i) = synthetic(i) + (-.05 + (.05+.05).*rand(1,1))*synthetic(i);
% end
% % %Generate synthetic data
% allExpParameters = xlsread('SyntheticData.xlsx','ExpParameters');
% % expParameters = allExpParameters;
% % synthetic = zeros(length(allExpParameters),1);
% % syntheticWithNoise = zeros(length(allExpParameters),1);
% % exitflag = zeros(length(allExpParameters),1);
% % syntheticPlasmaParams = zeros(plasmaUnknowns,length(allExpParameters));
% % % Calculate synthetic experiments
% % for i=1:length(allExpParameters)
% %     [synthetic(i),syntheticPlasmaParams(:,i),exitflag(i)]= GlobalSolver(real,i);
% % end
% % %Add noise to synthetic experiments
% % for i=1:length(allExpParameters)
% %     syntheticWithNoise(i) = synthetic(i) + normrnd(0,expError);
% % end
% % 
% 
% %Removes unrealistic etch rates
% % index = 1;
% % remove = [];
% % for i = 1:length(synthetic)
% %     if(synthetic(i)>.5e+9)
% %         remove(index) = i;
% %         index = index+1;
% %     end
% % end
% % synthetic(remove,:) = [];
% % allExpParameters(remove,:) = [];
% %load allExpParameters
% %Split data into training and test sets
% trainingDataIndex = randperm(size(synthetic,1)+1);
% trainingDataIndex = trainingDataIndex(1:size(synthetic,1)/6)-1;
% noTrainingCases = length(trainingDataIndex);
% noTestCases = size(synthetic,1) - noTrainingCases;
% trainingData = zeros(noTrainingCases,1);
% trainingExp = zeros(noTrainingCases,noExpParameters);
% testData = [];
% testExp = [];
% for i=1:noTrainingCases
%     trainingData(i) = synthetic(trainingDataIndex(i));
%     trainingExp(i,:) = allExpParameters(trainingDataIndex(i),:);
% end
% Lia = ismember(synthetic,trainingData);
% for i=1:length(Lia)
%     if ~(Lia(i))
%         testData = [testData synthetic(i)];
%         testExp = [testExp; allExpParameters(i,:)];
%     end
% end
% 
%Generate Test synthetic Data
x = linspace(1,10,5)
%Make initial guess for unknown parameters
current = zeros(noUnknowns,1);
for i = 1:3
          current = ProposeParameters(current,i);
end
trainingExp = x;
expParameters = trainingExp;
trainingData = zeros(length(trainingExp),1);
for i=1:length(trainingExp)
    trainingData(i) = testArr(real,i); 
end
%set experiments to training data

data = trainingData;

index = 1;
nn      = 100;       % Number of samples for examine the AC
N       = 1000;     % Number of samples (iterations)
burnin  = 1;      % Number of runs until the chain approaches stationarity
lag     = 1;        % Thinning or lag period: storing only every lag-th point
theta   = zeros(N,noUnknowns); 
acc = zeros(subBlocks,1);
[PosteriorCurrent] = Posterior(current,1);
count = 0;
AlphaSet = zeros(N,subBlocks);

% for i = 1:burnin    % First make the burn-in stage
%     [t] = MetropolisHastings(current);
% end
%Peform MH iterations
totalTime = tic;

for cycle = 1:N  % Cycle to the number of samples
    %for j = 1:lag 
    MHtime = tic;
    for j=1:3 % Cycle to make the thinning
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

for i =1:noUnknowns
    figure; hist(theta(:,i),100)
end
draws = N;
param = zeros(1,3);
pred = zeros(N,length(data));


for j = 1:draws
        param(1) = theta(round(rand(1)*N)+1,1);
        param(2) = theta(round(rand(1)*N)+1,2);
        param(3) = theta(round(rand(1)*N)+1,3);
    for i =1:length(expParameters)
        pred(j,i) = testArr(param,i);
    end
end

figure;
subplot(2,3,1);
hist(pred(:,1),100)
hold on
line([data(1) data(1)],[0 100], 'Color', 'r')
subplot(2,3,2),100;
hist(pred(:,2))
hold on
line([data(2) data(2)],[0 100], 'Color', 'r')
subplot(2,3,3);
hist(pred(:,3))
hold on
line([data(3) data(3)],[0 100], 'Color', 'r')
subplot(2,3,4);
hist(pred(:,4))
hold on
line([data(4) data(4)],[0 100], 'Color', 'r')
subplot(2,3,5);
hist(pred(:,5))
hold on
line([data(5) data(5)],[0 100], 'Color', 'r')

figure;
hist(pred(:,1),100)
hold on
line([1500 1500],[0 NaN], 'Color', 'r')
% figure; scatter(x,pred);hold on; scatter(x,data,'g')
% xlabel('Training Experiment')
% ylabel('Output')
% title('N = 500, Experiments with Known Error')
% %simulate training data
% simTrainingData = zeros(length(trainingData),1);
% simPlasmaParams = zeros(plasmaUnknowns,length(trainingData));
% for i=1:length(data)
%     [simTrainingData(i),simPlasmaParams(:,i),exitflag(i)] = GlobalSolver(mean(theta),i);
% end
% figure(1)
% x = linspace(1,length(trainingData),length(trainingData));
% scatter(x,simTrainingData,'r')
% hold on
% scatter(x,trainingData,'g')
% 
% %Simulate test data
% expParameters = testExp;
% simTestData = zeros(length(testData),1);
% simTestPlasmaParams = zeros(plasmaUnknowns,length(testData));
% for i=1:length(testExp)
%     [simTestData(i),simTestPlasmaParams(:,i),extiflag] = GlobalSolver(mean(theta),i);
% end
% 
% figure(2)
% x = linspace(1,length(testData),length(testData));
% scatter(x,simTestData,'r')
% hold on
% scatter(x,testData,'g')


