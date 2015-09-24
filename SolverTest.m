clear
%Test function for global plasma system
%All SI units (m)

global massIon %ion mass
global kb %boltzmann constant
global q %electron charge
global R %radius of chamber
global L %length of chamber
global T %temperature of system (treat constant)
global V %volume
global A %area
global K %this is equal to k7 in Efremov study
global plasmaUnknowns;
global k9

%Define Constants
massIon = 35.45*(1.660568e-27); %kg
kb = 1.381e-23; %J/K
q = 1.6022e-19; %C
R = 0.15; %m
L = 0.14; %m
T = 303;
V = pi*R^2*L;
A = 2*pi*R*L + 2*pi*R^2;
K = 5e-14; %(m^3/2)
plasmaUnknowns = 14;

%Calculate k9
diffusionLength = sqrt(1/(2.405/R)^2+(pi/L)^2)
Df = diffusionLength/3*sqrt(8*kb*T/(pi*massIon)); %need to calculate effective diffusion coefficient
recombinationProbability = (0.2+10e-3+0.05)/3;
k9 = recombinationProbability*Df/(diffusionLength^2);

%Define means for Bayesian parameters
k7_nond = 5e-8;
k1 = 3e-10/k7_nond;
k2 = 2.1e-12/k7_nond;
k3 = 1.5e-10/k7_nond;
k4 = 3e-9/k7_nond;
k5 = 1e-10/k7_nond;
k6 = 2e-11/k7_nond;
k8 = k2;
mean = [k1 k2 k3  k4 k5 k6 k8 2 6 3 2 4 1 2];
Act = mean(1:7);
B = mean(8:14);

%Define intitial guess
parameter(1).name = 'nCl2';
parameter(2).name = 'nCl';
parameter(3).name = 'nAr';
parameter(4).name = 'ne';
parameter(5).name = 'nCl_neg';
parameter(6).name = 'nCl2_pos';
parameter(7).name = 'nCl_pos'; 
parameter(8).name = 'nAr_pos';
parameter(9).name ='Te';
parameter(10).name = 'Beta';
parameter(11).name = 'BetaS';
parameter(12).name = 'gamma_T';
parameter(13).name = 'dc';
parameter(14).name = 'v';
x0 = [1 1 1 1 1 1 1  1 3 1 1 20 A 10e+6];

%Solve initial system
opts = optimset('MaxFunEvals',10e+8,'MaxIter',10e+8,'display','iter');
%x1 = fsolve(@GlobalPlasmaSystem,x0,opts);
opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','off','TolCon',10e-10);
l = zeros(plasmaUnknowns,1);
u = [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf];
x1 = fmincon(@(x)0,x0,[],[],[],[],l,u,@(x)fminconstr(x,Act,B,1),opts);
x2 = MakeDimensional(x1,1);
