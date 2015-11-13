function [etchRate,xFinal,exitflag] = GlobalSolver(current,expNo)
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
global sigma
global expParameters;
P = expParameters(expNo,1);
n0 = P/(kb*T);
Te = expParameters(expNo,7); %for synthetic experiment generation
kGuess = 4e-16;

%Grab parameter guesses
Act = current(1:7);
B = current(8:14);

%Normalizations
lambdaMax = 1/(sigma*n0);
maxTe = 25;
vMax = sqrt(q*maxTe/massIon);
DMax = lambdaMax*vMax;
kNorm = 10e-20;
nNorm = n0;
x0 = [
n0/nNorm;
n0/nNorm;
n0/nNorm;
n0/nNorm;
n0/nNorm;
n0/nNorm;
n0/nNorm;
n0/nNorm;
1
1
1
1
.01
10e+3/vMax;
1
1
1
1
kGuess/kNorm;
kGuess/kNorm;
kGuess/kNorm;
kGuess/kNorm;
kGuess/kNorm;
kGuess/kNorm;
kGuess/kNorm];

%Reasonable unormalized values for reference
Beta = 1.6;
BetaS = 7.44e-1;
gamma_T = 28;
v = 70610;
nCl2 = 7.57e+19;
nAr = 1.21e+20;
ne = 8e+19;
nCl_neg = 1.28e+17;
nCl2_pos = 1.3e+17;
nCl_pos = 9e+17;
nAr_pos = 3e+17;
dc = 0.0729;
lambda = 1/(sigma*n0);




l = zeros(plasmaUnknowns,1);
u = Inf(plasmaUnknowns,1);
l(10) = 1;
l(12) = 1;
u(1) = 1;
u(2) = 1;%nCl;
u(3) = 1;
u(4) = 1/1000;
u(5) = 1/1000;
u(6) = 1/1000;
u(7) = 1/1000;
u(8) = 1/1000;
u(9) = 5;
u(10) = 4;
u(14) = 1;
u(15) = 1;
u(16) = 1;


opts = optimoptions(@fmincon,'TolCon',10e-3,'MaxFunEvals',10e+4,'Maxiter',10e+4,'TolX',10e-7);
problem = createOptimProblem('fmincon','objective',@(x)0,'nonlcon',@(x)fminconstr(x,Act,B,expNo),'x0',x0,'lb',l,'ub',u,'options',opts);
trials = 20;
startPoints = zeros(trials,plasmaUnknowns);
for i=1:trials
    for n=1:8
        startPoints(i,n) = 1/i+10e-4;%unifrnd(10e-3,1000);
    end
    for n=9:18
        startPoints(i,n) = 20/i+10e-3;%unifrnd(10e-3,20);
    end
    startPoints(i,10) = 20/i+1;%unifrnd(1,20);
    startPoints(i,12) = 20/i+1; %unifrnd(1,20);
    for n=19:25
        startPoints(i,n) = 1000/i+10e-5;%unifrnd(10e-5,10);
    end
end
tpoints = CustomStartPointSet(startPoints);
MS = MultiStart;
[x,fval,exitflag]  = run(MS,problem,tpoints);
xFinal = Rescale(x,expNo);
if exitflag==2
etchRate = CalcEtchRate(xFinal,expNo);
else
   etchRate = 10e+8; %picking unrealistic values so that these parameters are not chosen
end

%F = [etchRate,xFinal];
