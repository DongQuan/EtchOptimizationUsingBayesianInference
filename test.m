global L
L = .14;
current = ones(15,1);
current = current*10e-5;
Act = current(1:7)
B = current(8:14);
x = ones(25,1);
%check = GlobalPlasmaSystemWithDimensionsAndK(x,Act,B,expNo)
%s = GlobalSolver(current,1);
m = Likelihood(current);