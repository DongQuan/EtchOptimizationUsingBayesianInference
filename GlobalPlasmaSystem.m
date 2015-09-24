function F = GlobalPlasmaSystem(x,A2,B,i)
global ExpParameters;
global m;
global kb;
global q_e
global TeSet
global R
global L
global T
i=1;
%Constants
A = 2*pi*R*L + 2*pi*R^2;
V_nond = pi*R^2*L;
kb = 1.381*10^-23;
m = 35.45*1.660468e-27; 
k7_nond = 5e-8;
%PrevX = MakeDimensional(x,i);
%Te = PrevX(10);
%Experimental Parameters
%Te = TeSet(1);
k7_nond = 5e-8;
k1 = 3e-10/k7_nond;
k2 = 2.1e-12/k7_nond;
k3 = 1.5e-10/k7_nond;
k4 = 3e-9/k7_nond;
k5 = 1e-10/k7_nond;
k6 = 2e-11/k7_nond;
k8 = k2;
% k = CalculateRates(A2,B,Te);
% assignin('base', 'k', k*k7_nond);
% k1 = k(1);
% k2 = k(2);
% k3 = k(3);
% k4 = k(4);
% k5 = k(5);
% k6 = k(6);
% k8 = k(7);
%P, T, Flowrate, PIcp, Prf
P0_nond = ExpParameters(i,1);
T0_nond = T;
FlowRate = ExpParameters(i,2);
FlowRate = FlowRate/60; %divide by 60 for per second
tau_nond = V_nond/FlowRate;


n0_nond = P0_nond/(kb*T0_nond)/(100^3);%convert to cm^3

k9 = 2*10^3/(n0_nond*k7_nond);

%Non-dimensionalizations
n0 = 1;
k7 = 1;
tau = tau_nond*k7_nond*n0_nond;
Picp = ExpParameters(i,3);
Prf = ExpParameters(i,4);
W_nond = (Prf+Picp)*(100^2);
delta = ExpParameters(i,5);
dc = 7.3;
T_dim = m/q_e*(k7_nond*n0_nond*x(11))^2;
Te = 2;
Te = Te/T_dim; 

%Functions
F=[
2*(k4+ k3)*(x(4)) *(x(1))-(k9+1/tau)*(x(2)); %1, (2*k4+k3)*ne*nCl2-(x(2)^2)*(k9+1/tau_dim)Te
(x(1))-(1-delta)*(n0)+ 0.5*(x(2)); %2, nCl2-n0*(1-delta)*T0/T+0.5*nCl 
(x(3))- delta*n0;%3,%nAr- n0*delta*T0/T %checked through here
(k2+k3)*(x(1))*((x(4)))-k7*((x(7))+(x(8))+(x(6)))*(x(5))-Te*(x(5))*x(4); %4, (k3+k3)*nCl2*ne-k7*(nCl_pos+nAr_pos+nCl2_pos)*nCl_neg-k8*nCl_neg*ne %nCl2+
k1*(x(1))*((x(4)))-((x(12)) + k7*((x(5))))*(x(6)); % 5, k1*nCl2*(x(4)^2)-(v/dc + k7*(x(5)^2))*(x(6)^2)x(1
(k5*((x(2)))+k2*(x(1)))*((x(4)))-((x(12)) + k7*((x(5))))*(x(7));%6 k5, 23
k6*(x(3))*((x(4)))-((x(12))+ k7*((x(5))))*(x(8)); %7 x(7, k6*nAr*ne-(v/dc  + k7*nCl_neg)*nAr_pos %nAr+
((x(5)))+(x(4))- (x(6))-(x(7))-(x(8)); %8 quasineutrality
((x(5)))/((x(4))) - (x(10))*exp((1+(x(10)))*(x(13)-1)/(2*(1+(x(10))*x(13))));%9, Beta - BetaS*exp((1+BetaS)*(gamma_T-1)/(2*(1+BetaS*gamma_T)))
(x(12)) - sqrt((Te))*sqrt(1+(x(10))/(1+(x(10))*(x(13))));%10, v- sqrt((1+BetaS)/(1+BetaS*gamma_T))*sqrt(q_e*Te/m)/(k7*n0*R*L/(R*hl+L*hr))
2*k1*((x(1))+k5*((x(2)))+k6*(x(3)))+(k2-k3)*(x(1))-(1+((x(5)))/((x(4))))*(2*((x(12))+k7*((x(5)))/((x(4)))*((x(4)))))+Te*((x(4)))*((x(5)))/((x(4)));%11;%17, 2*(k1*nCl2+k5*nCl+k6*nAr)+(k3-k3)*nCl2-(1+Beta)*(2*v/dc+k7*Beta*ne)+x(10)*ne*Beta
W_nond/(k7_nond^3*n0_nond^4*((x(11))^3)*A*m)-(x(4))*(x(12))*(7*(Te));%Check this eqn
x(11)-(0.5*R*L/(R*( (1+2*x(10)/x(13))/(1+x(10)))+L*(1+3*x(10)/x(13))/(1+x(10))));];%12

 
% % 

end