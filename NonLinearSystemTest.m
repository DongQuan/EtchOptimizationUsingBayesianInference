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
global epMatrix

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

kNorm = 10e+10;
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
a1 = 5;
a2 = 3;
a3 = 2;
real = [a1 a2 a3];
center = [.8 5 1];
sd = [.5 2 .25]
%sd = ones(15,1)*.25;
for i = 1:length(center)
    m = center(i);
    v = sd(i)^2;
    proposalSD(i) = sqrt(log(v/(m^2)+1));
    proposalCenter(i) = log((m^2)/sqrt(v+m^2));
end

noExp = 3
epMatrix = rand(3,3);
%Make initial guess for unknown parameters
current = zeros(1,noUnknowns);
for i = 1:2
          current = ProposeParameters(current,i);
end
trainingExp = epMatrix;
expParameters = trainingExp;
trainingData = zeros(length(trainingExp),1)
for i=1:noExp
    trainingData(i) = testMatrix(real,i); 
%set experiments to training data
end
data = trainingData;

index = 1;
nn      = 100;       % Number of samples for examine the AC
N       = 500;     % Number of samples (iterations)
burnin  = 1;      % Number of runs until the chain approaches stationarity
lag     = 1;        % Thinning or lag period: storing only every lag-th point
theta   = zeros(N*2,noUnknowns); 
acc = zeros(subBlocks,1);
[PosteriorCurrent] = Posterior(current,1);
count = 0;
AlphaSet = zeros(N,subBlocks);

% for i = 1:burnin    % First make the burn-in stage
%     [t] = MetropolisHastings(current);
% end
%Peform MH iterations
totalTime = tic;
current(3) = .1;
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
% 
% for i =1:noUnknowns
%     figure; hist(theta(:,i),100)
% end
% draws = 1000;
% param = zeros(1,3);
% param(3) = 3;
% 
% for i =1:length(expParameters)
%     pred(i) = testArr(mean(theta),i);
% end
% for i =1:length(expParameters)
%     for j = 1:draws
%         param(1) = theta(round(rand(1)*N),1);
%         param(2) = theta(round(rand(1)*N),2);
%         pred(j,i) = testArr(param,i);
%     end
% end
% 
% figure;
% subplot(2,3,1);
% hist(pred(:,1),100)
% subplot(2,3,2),100;
% hist(pred(:,2))
% subplot(2,3,3);
% hist(pred(:,3))
% subplot(2,3,4);
% hist(pred(:,4))
% subplot(2,3,5);
% hist(pred(:,5))
% 
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
% 
% 
