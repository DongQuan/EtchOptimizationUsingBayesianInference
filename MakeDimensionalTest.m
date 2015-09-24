function f = MakeDimensionalTest(x3,i)
global PlasmaUnknowns
global kb;
global q_e;
global m;
global ExpParameters;
%Constants



%Experimental Parameters
P0_nond = ExpParameters(i,1);
T0_nond = ExpParameters(i,2);
n0_nond = P0_nond/(kb*T0_nond)/(100^3);%convert to cm^3

dc = 7.3;
%Non-dimensionalizations
k7_nond = 5e-8;
Te = 2.3;
%Test Val
x4=zeros(PlasmaUnknowns,1);
T_dim = m/q_e*(k7_nond*n0_nond*dc)^2;
assignin('base', 'T_dim', T_dim);
for i=1:8
    x4(i) = x3(i)*n0_nond;
end
x4(9) = x3(9)*T0_nond;
x4(10) = Te;
x4(12) = x3(12);
x4(13) = x3(13);
x4(14) = dc;
x4(15) = x3(15)*(k7_nond*n0_nond*dc);
x4(16) = x3(16);
x4(17) = x3(17);
x4(11) = Te/x4(16);
f = x4;