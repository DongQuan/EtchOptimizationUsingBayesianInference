function f = MakeDimensional(x1,i)
global plasmaUnknowns
global kb
global q
global massIon
global K
global T
%Experimental Parameters
P = 2; %Pa
n0 = P/(kb*T)
dcConstant = 10; %This constant is a placeholder for when we solve for hl and hr
x2=zeros(plasmaUnknowns,1);

for i=1:8
    x2(i) = x1(i)*n0;
end
x2(10) = x1(10);
x2(11) = x1(11);
x2(12) = x1(12);
x2(13) = x1(13)*0.5*R*L/(R*hl+l*hr);
T_dim = massIon/q*(K*n0*x2(13))^2;
assignin('base', 'T_dim', T_dim);
x2(9) = x1(9)*T_dim*dcConstant;
x2(14) = x1(13)*(K*n0*x2(11))

f = x2;