%function F = GlobalSolver(current,expNo)
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
T = 303;
n0 = P/(kb*T);




opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','off','TolCon',10e-10);
l = zeros(plasmaUnknowns,1);
l(9) = 2;
% l(14) = 50;
%Define bounds for rates
k1 = 3e-16;
k2 = 2.1e-18;
k3 = 1.5e-16;
k4 = 3e-15;
k5 = 1e-16;
k6 = 2e-17;
k8 = k2;
kGuess = 4e-16;
current = [k1 k2 k3  k4 k5 k6 k8 0 0 0 0 0 0 0 0];
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

Te = 2.3;
Beta = 1.6;
BetaS = 7.44e-1;
gamma_T = 28;
v = 70610;
nCl2 = 7.57e+19;
nCl = 9e+19;
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
% l(9) = Te;
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
% u(9) = 2.3;
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
u = Inf(plasmaUnknowns,1);
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
opts = optimoptions(@fmincon,'TolCon',10e-10,'MaxFunEvals',10e+4,'Maxiter',10e+4);
problem = createOptimProblem('fmincon','objective',@(x)0,'nonlcon',@(x)fminconstr(x,Act,B),'x0',x0,'lb',l,'ub',u,'options',opts);
trials = 10; %4:07
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
%[x1,f1] = fmincon(problem);
tpoints = CustomStartPointSet(startPoints);
MS = MultiStart;
%gs = GlobalSearch(ms,'NumTrialPoints',1000,'TolX',10e-20);
[x, f] = run(MS,problem,tpoints);
xFinal = Rescale(x);
%[xg,gf,exitflag,output,solutions]=run(gs,problem)

%xg = fmincon(@(x)0,x0,[],[],[],[],l,u,@(x)fminconstr(x,Act,B),opts);

PredER = CalcEtchRate(xFinal,31);


%F = xg;
