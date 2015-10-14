function F = GlobalSolver(current,expNo)
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

%Grab parameter guesses
Act = current(1:7);
B = current(8:14);



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
kGuess = 10e-16;
%x = @(x)fminconstr(x,Act,B,expNo);
x0 = [
10E+17
10E+17
10E+17
10E+17
10E+14
10E+17
10E+17
10E+17
1
1
1
1
.01
10e+3
1
1
1
1
k1
k2
k3
k4
k5
k6
k8];

% x0 = [
% n0
% n0
% n0
% n0
% n0
% n0
% n0
% n0
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
% kGuess
% kGuess
% kGuess
% kGuess
% kGuess
% kGuess
% kGuess];

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


%x0 = ones(25,1);

%u = [n0 n0 n0 n0/1000 n0/100 n0/100 n0/100 n0/100 2 Inf Inf Inf Inf 2000 Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf];
l = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 k1  k2 k3 k4 k5 k6 k8];
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
u = [Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf Inf  Inf Inf Inf Inf Inf  k1  k2 k3 k4 k5 k6 k8];
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
%fun = @(x)GlobalPlasmaSystemWithDimensionsAndK(x,Act,B);
%fun = @(x)0,x0,[],[],[],[],l,[],@fminconstr;
%problem = createOptimProblem('fmincon','objective',@(x)0,x0,[],[],[],[],l,[],@fminconstr,'x0',x0);
%problem = createOptimProblem('fmincon','objective',fun,'x0',x0,'lb',l,'ub', u);
%stpoints = RandomStartPointSet('NumStartPoints',10);
%ms = MultiStart;
%[xmulti, fval, exitflag, output]=run(ms,problem,stpoints)
opts = optimoptions(@fmincon,'TolCon',10e-3);
problem = createOptimProblem('fmincon','objective',@(x)0,'nonlcon',@(x)fminconstr(x,Act,B),'x0',x0,'lb',l,'ub',u,'options',opts);
%[x1,f1] = fmincon(problem);
gs = GlobalSearch('NumTrialPoints',10000,'TolX',10e-14);
[xg,gf,exitflag,output,solutions]=run(gs,problem)
%[xmulti, fval] = gamultiobj(fun,25,[],[],[],[],l,u)
%PredER = CalcEtchRate(xmulti,expNo);


F = xg;
