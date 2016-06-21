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
global priorSD
global priorUB
global priorLB
global lb
global ub
global proposalLB
global proposalUB
global priorCenter;
global kNorm

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
noUnknowns = 14; %7 k's (A and B coeff) plus standard error
noExpParameters = 9;
kNorm = 10e+20;

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


proposalLB = zeros(1,noUnknowns);
proposalUB = [1000 1000 1000 1000 1000 1000 1000 20 20 20 20 20 20 10];
priorLB = zeros(1,noUnknowns);
priorUB = [500 500 2000 1000 300  1000 500 20 20 20 20 20 20 10];



synthetic = [10.57496248	4.558329124	4.293246598	18.13634377	10.40168263	8.534352014	5.352673677	9.60319921	13.82963681	18.49105692	4.036717427	14.80609507	9.141620243	8.98011037	14.60610744	16.05577287	13.91058772	12.13618363	0.988809673	0	20.15489026	30.53190617	2.465218094	0	9.296565152	0.00066889	3.92E-08	1.96058308	15.07131699]
allExpParameters = FactorialDesign();
% %load allExpParameters
% %Split data into training and test sets
trainingDataIndex = randperm(length(synthetic)+1);
trainingDataIndex = trainingDataIndex(1:round(length(synthetic)/7))-1;
noTrainingCases = length(trainingDataIndex);
noTestCases = size(synthetic,1) - noTrainingCases;
trainingData = zeros(noTrainingCases,1);
trainingExp = zeros(noTrainingCases,noExpParameters);
testData = [];
testExp = [];
for i=1:noTrainingCases
    trainingData(i) = synthetic(trainingDataIndex(i));
    trainingExp(i,:) = allExpParameters(trainingDataIndex(i),:);
end
Lia = ismember(synthetic,trainingData);
for i=1:length(Lia)
    if ~(Lia(i))
        testData = [testData synthetic(i)];
        testExp = [testExp; allExpParameters(i,:)];
    end
end


%set experiments to training data
data = trainingData;
expParameters = trainingExp;


index = 1;
nn      = 100;       % Number of samples for examine the AC
N       = 1;     % Number of samples (iterations)
burnin  = 500;      % Number of runs until the chain approaches stationarity
lag     = 1;        % Thinning or lag period: storing only every lag-th point
theta   = zeros(N*noUnknowns,noUnknowns); 
acc = zeros(noUnknowns,1);

count = 0;
AlphaSet = zeros(N,noUnknowns);


%Peform MH iterations
totalTime = tic;
%Make initial guess for unknown parameters
current = zeros(noUnknowns,1);
for i = 1:7
          current = ProposalFunction(theta,current,i);
end
[PosteriorCurrent] = Posterior(current,1);
ind = 0;

    MHtime = tic;
    for cycle=1:N
       for m=1:7 % Cycle to make the thinning
            SCtime = tic;
            [alpha,t, a,prob, PosteriorCatch] = MetropolisHastings(theta,current,PosteriorCurrent,m);
            SCelapsed = toc(SCtime);
            theta((cycle-1)*7+m,:) = t;        % Samples accepted
            AlphaSet(cycle,m) = alpha;
            current = t;
            PosteriorCurrent = PosteriorCatch;
            acc(m) = acc(m) + a;  % Accepted ?
        end
       %plotter(theta(ind*N+1:(ind+1)*N,:),N,trainingExp,testExp,trainingData,real)
    end
    count = count+1;
    MHelapsed = toc(MHtime);
%end
totalTimElapsed = toc(totalTime);
accrate = acc/N;     % Acceptance rate,. 

%plotter(theta,N*noUnknowns,trainingExp,testExp,trainingData,testData)



figure; 
xlabel('Training Experiment')
ylabel('Output')
title('N = 100, Experiments with Known Error')
%simulate training data
simTrainingData = zeros(length(trainingData),1);
simPlasmaParams = zeros(plasmaUnknowns,length(trainingData));
for i=1:length(data)
    [simTrainingData(i),simPlasmaParams(:,i),exitflag(i)] = GlobalSolver(mean(theta),i);
end
figure(1)
x = linspace(1,length(trainingData),length(trainingData));
scatter(x,simTrainingData,'r')
hold on
scatter(x,trainingData,'g')

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


