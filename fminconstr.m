function [c,ceq] = fminconstr(x,Act,B,expNo)
global massIon
global expParameters
global kb
global T
global q

P = expParameters(expNo,1); %Pa
n0 = P/(kb*T);%initial neutral density
maxTe = 25;
vMax = 1/(sqrt(q*maxTe/massIon));
nNorm = 1/n0;

%Surface kinetics parameters
RxnProb = .25;
eos = 40; %eV
eod = 10; % eV
A1 = .05;
B1 = A1;
density_MgO = 3.56e+3; %kg/m^3
M = 40.3;
Na = 6.02e+23;
c=[6e+7*M/(Na*density_MgO)*(1-RxnProb*((x(2)/nNorm)*(x(14)/vMax))/(RxnProb*((x(2)/nNorm)*(x(14)/vMax)) + (B1*(sqrt((6*x(9)))-sqrt(eod)))*((x(14)/vMax)*((x(6)/nNorm)+(x(7)/nNorm)+(x(8)/nNorm)))))*(RxnProb*((x(2)/nNorm)*(x(14)/vMax))+(A1*(sqrt((6*x(9)))-sqrt(eos)))*((x(14)/vMax)*(x(8)/nNorm)))-200];
ceq = GlobalPlasmaSystemWithDimensionsAndK(x,Act,B,expNo); % the fsolve objective is fmincon constraints