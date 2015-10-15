
k1 = 3e-10;
k2 = 2.1e-13;
k3 = 1.5e-10;
k4 = 3e-9;
k5 = 1e-10;
k6 = 2e-11;
k8 = k2;
noise = 10;
expParameters = FactorialDesign();
current = [k1 k2 k3  k4 k5 k6 k8 0 0 0 0 0 0 0 0];
%Grab parameter guesses
Act = current(1:7);
B = current(8:14);

x0 = ones(8,1);
l = zeros(8,1);
u = [Inf Inf Inf Inf Inf Inf Inf Inf];
opts = optimoptions(@fmincon,'TolCon',10e-12);
problem = createOptimProblem('fmincon','objective',@(x)0,'nonlcon',@(x)fminconstr(x,Act,B),'x0',x0,'lb',l,'ub',u,'options',opts);
%[x1,f1] = fmincon(problem);
gs = GlobalSearch('NumTrialPoints',10000,'TolX',10e-14);
[xg,gf,exitflag,output,solutions]=run(gs,problem)