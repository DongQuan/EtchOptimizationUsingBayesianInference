clear

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
T = 600;
V = pi*R^2*L;
A = 2*pi*R*L + 2*pi*R^2;
K = 5e-14; %(m^3/s)
plasmaUnknowns = 25;
sigma = 16.8e-24;
noUnknowns = 15; %7 k's (A and B coeff) plus standard error
%Calculate k9
diffusionLength = sqrt(1/(2.405/R)^2+(pi/L)^2);
Df = diffusionLength/3*sqrt(8*kb*T/(pi*massIon)); %need to calculate effective diffusion coefficient
recombinationProbability = (0.2+10e-3+0.05)/3;
k9 = recombinationProbability*Df/(diffusionLength^2);

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
current = [k1 k2 k3  k4 k5 k6 k8 2 2 2 2 2 2 2];
%current = [k1 k2 k3  k4 k5 k6 k8 2 5 6 10 18 5 8];
xmulti=zeros(plasmaUnknowns,length(expParameters));
predER = zeros(1,10);
for expNo=1:3%:length(expParameters)
        start = tic;
        [predER(expNo),xmulti(:,expNo),flag] = GlobalSolver(current,expNo);
        elapsed = toc(start);
        %PredER(expNo) = CalcEtchRate(xmulti(:,expNo),expNo);
end
%SyntheticDataWithNoise = SyntheticData + noise*randn(length(SyntheticData));


