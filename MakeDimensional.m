function f = MakeDimensional(x3,i)
global PlasmaUnknowns
global kb;
global q_e;
global m;
global ExpParameters;
global R
global L

%Experimental Parameters
P0_nond = ExpParameters(i,1);
T0_nond = ExpParameters(i,2);
n0_nond = P0_nond/(kb*T0_nond)/(100^3);%convert to cm^3


%Non-dimensionalizations
k7_nond = 5e-8;
%Test Val
x4=zeros(PlasmaUnknowns,1);
T_dim = m/q_e*(k7_nond*n0_nond*x3(11))^2;
assignin('base', 'T_dim', T_dim);
for i=1:8
    x4(i) = x3(i)*n0_nond;
end
x4(9) = x3(9)*T_dim;
x4(10) = x3(10);
x4(11) = x3(11);
x4(12) = x3(12)*(k7_nond*n0_nond*x4(11));
x4(13) = x3(13);

f = x4;