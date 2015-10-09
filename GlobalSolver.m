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


%Grab parameter guesses
Act = current(1:7)
B = current(8:14)




opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','off','TolCon',10e-10);
l = zeros(plasmaUnknowns,1);
l(9) = .1;
l(14) = 50;
%Define bounds for rates
k1 = 3e-16;
k2 = 2.1e-18;
k3 = 1.5e-16;
k4 = 3e-15;
k5 = 1e-16;
k6 = 2e-17;
k8 = k2;

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



P = expParameters(expNo,1);
n0 = P/(kb*T);
u = [n0 n0 n0 n0/1000 n0/100 n0/100 n0/100 n0/100 15 Inf Inf Inf Inf 2000 Inf Inf Inf Inf 10e-14 10e-14 10e-14 10e-14 10e-14 10e-14 10e-14];
fun = @(x)GlobalPlasmaSystemWithDimensionsAndK(x,Act,B);
problem = createOptimProblem('lsqnonlin','objective',fun,'x0',x0,'lb',l,'ub', u);
stpoints = RandomStartPointSet('NumStartPoints',10,'ArtificialBound',200);
ms = MultiStart;
[xmulti]=run(ms,problem,stpoints)
PredER = CalcEtchRate(xmulti,expNo);


F = xmulti;
