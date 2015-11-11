function [etchRate,xFinal] = GlobalSolver(current,expNo)
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
Te = expParameters(expNo,7);
kGuess = 4e-16;
etchGuess = 50;

%Grab parameter guesses
Act = current(1:7);
B = current(8:14);

% x0 = [
% 10E+17
% 10E+17
% 10E+17
% 10E+17
% 10E+14
% 10E+17
% 10E+17
% 10E+17
% 1
% 1
% 1
% 1
% .01
% 10e+3
% 1
% 1
% 1
% 1
% k1
% k2
% k3
% k4
% k5
% k6
% k8];
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


%x0 = ones(plasmaUnknowns,1);


%u = [n0 n0 n0 n0/1000 n0/100 n0/100 n0/100 n0/100 2 Inf Inf Inf Inf 2000 Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf];

% l(1) = nCl2;
% l(2) = nCl;
% l(3) = nAr;
% l(4) = ne;
% l(5) = nCl_neg;
% l(6) = nCl2_pos;
% l(7) = nCl_pos;
% l(8) = nAr_pos;
%l(9) = Te;
% l(10) = Beta;
% l(11) = BetaS;
% l(12) = gamma_T;
% l(13) = dc;
% l(14) = v;
% l(16) = lambda;

% u(1) = nCl2;
% u(2) = nCl;
% u(3) = nAr;
% u(4) = ne;
% u(5) = nCl_neg;
% u(6) = nCl2_pos;
% u(7) = nCl_pos;
% u(8) = nAr_pos;
%u(9) = Te;
% u(10) = Beta;
% u(11) = BetaS;
% u(12) = gamma_T;
% u(13) = dc;
% u(14) = v;
% u(16) = lambda;


% l(19) = k1/kNorm;
% l(20) = k2/kNorm;
% l(21) = k3/kNorm;
% l(22) = k4/kNorm;
% l(23) = k5/kNorm;
% l(24) = k6/kNorm;
% l(25) = k8/kNorm;
l = zeros(plasmaUnknowns,1);
l(10) = 1;
l(12) = 1;
%l(2) = nCl;
u = Inf(plasmaUnknowns,1);
%l(9) = Te;
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

% l(19) = k1;
% l(20) = k2;
% l(21) = k3;
% l(22) = k4;
% l(23) = k5;
% l(24) = k6;
% l(25) = k8;
% u(19) = k1;
% u(20) = k2;
% u(21) = k3;
% u(22) = k4;
% u(23) = k5;
% u(24) = k6;
% u(25) = k8;
% u(19) = k1/kNorm;
% u(20) = k2/kNorm;
% u(21) = k3/kNorm;
% u(22) = k4/kNorm;
% u(23) = k5/kNorm;
% u(24) = k6/kNorm;
% u(25) = k8/kNorm;
%u(26) = 200;
opts = optimoptions(@fmincon,'TolCon',10e-3,'MaxFunEvals',10e+4,'Maxiter',10e+4,'TolX',10e-6);
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
    %startPoints(i,26) = 50/i+1; %unifrnd(1,20);
end
%[x1,f1] = fmincon(problem);
tpoints = CustomStartPointSet(startPoints);
MS = MultiStart;
%gs = GlobalSearch(ms,'NumTrialPoints',1000,'TolX',10e-20);
[x, f] = run(MS,problem,tpoints);

%[xg,gf,exitflag,output,solutions]=run(gs,problem)
%x = fmincon(@(x)0,x0,[],[],[],[],l,u,@(x)fminconstr(x,Act,B,expNo),opts);
xFinal = Rescale(x,expNo);
etchRate = CalcEtchRate(xFinal,expNo);

%F = [etchRate,xFinal];
