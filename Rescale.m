function F = Rescale(x1,expNo)
global sigma
global massIon
global expParameters
global kb
global T
global q
x2 = x1;
P = expParameters(expNo,1); %Pa
n0 = P/(kb*T);%initial neutral density
lambdaMax = 1/(sigma*n0);
maxTe = 25;
vMax = sqrt(q*maxTe/massIon);
DMax = lambdaMax*vMax;
kNorm = 10e-18;
nNorm = n0;
for i=1:8
    x2(i) = x1(i)*nNorm;
end
for i = 19:25
    x2(i) = kNorm*x1(i);
end
x2(14) = x1(14)*vMax;
x2(15) = x1(15)*DMax;
x2(16) = x1(16)*lambdaMax;
F = x2;