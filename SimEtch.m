function [F ,PlasmaParams] = SimEtch(current,i,expected)
global PlasmaUnknowns;
global ExpParameters;
global PlasmaSet;
P_nond = ExpParameters(i,1);
%x0 = ones(PlasmaUnknowns,1);
parameter(1).name = 'nCl2';
parameter(2).name = 'nCl';
parameter(3).name = 'nAr';
parameter(4).name = 'ne';
parameter(5).name = 'nCl_neg';
parameter(6).name = 'nCl2_pos';
parameter(7).name = 'nCl_pos'; 
parameter(8).name = 'nAr_pos';
parameter(9).name ='Te';
parameter(10).name = 'BetaS';
parameter(11).name = 'dc';
parameter(12).name = 'v';
parameter(13).name = 'gamma_T';

x0 = MakeDimensionless(expected,i);
A = current(1:7);
B = current(8:14);
% %Parameters that can be measured are: T
% %Parameters with bounds: Te, Ti, BiasFactor, 

opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','off','TolCon',10e-10);
l = zeros(PlasmaUnknowns,1);
diff_scale1 = 30/P_nond;
diff_scale2 = 30/(P_nond^2);
%l(14) = min(diff_scale1,diff_scale2);
%l(16) = 10/P_nond;
%l(10) = 9e-10;
u = [Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf,Inf];
%u(10) =10e-5;
%u(14) = diff_scale2;
x3 = fmincon(@(x)0,x0,[],[],[],[],l,u,@(x)fminconstr(x,A,B,i),opts);
opts = optimset('MaxFunEvals',10e+8,'MaxIter',10e+8);
%x3=fsolve(@GlobalPlasmaSystem,x0,opts);
assignin('base', 'fminsoln', x3);
x4 = MakeDimensional(x3,i);
PlasmaParams = x4;
PlasmaSet(:,i) = x4;
assignin('base', 'parameter', x4);
EtchRate = CalcEtchRate(i,x4);
 F = EtchRate;
end