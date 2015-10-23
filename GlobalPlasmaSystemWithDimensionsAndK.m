function F = GlobalPlasmaSystemWithDimensionsAndK(x,Act,B)
global massIon %ion mass
global kb %boltzmann constant
global q %electron charge
global T %temperature of system (treat constant)
global V %volume
global K %this is equal to k7 in Efremov study
global k9
global A
global L
global sigma
global R
global expParameters
%Constants
%PrevX = MakeDimensional(x,i);
expNo = 31;
Te = expParameters(expNo,7);


% k1 = 3e-16;
% k2 = 2.1e-18;
% k3 = 1.5e-16;
% k4 = 3e-15;
% k5 = 1e-16;
% k6 = 2e-17;
% k8 = k2;


%Define experimental parameters (with units)
P = expParameters(expNo,1); %Pa
flowRate = expParameters(expNo,2)/(100^3)/60; %m^3/min converted from sccm
tau = V/flowRate; %residence time (s)
delta = expParameters(expNo,5);
n0 = P/(kb*T);%initial neutral density
Picp = expParameters(expNo,3);
Prf = expParameters(expNo,4);
W = Picp + Prf;
TeV =  8.621738*10e-5 *T;
% parameter(1).name = 'nCl2';
% parameter(2).name = 'nCl';
% parameter(3).name = 'nAr';
% parameter(4).name = 'ne';
% parameter(5).name = 'nCl_neg';
% parameter(6).name = 'nCl2_pos';
% parameter(7).name = 'nCl_pos'; 
% parameter(8).name = 'nAr_pos';
% parameter(9).name ='Te';
% parameter(10).name = 'Beta';
% parameter(11).name = 'BetaS';
% parameter(12).name = 'gamma_T';
% parameter(13).name = 'dc';
% parameter(14).name = 'v';
% parameter(15).name = 'D';
% parameter(16).name = 'lambda';
% parameter(17).name = 'hl';
% parameter(18).name = 'hr';

%Normalization: Function variables are normalized according to their maxed
%values. Parameters 1-8 are normalized by n0. v (parameter 14) is normalized by
%sqrt(e*maxTe/m). Lambda is normalized by 1/(sigma*n0). D (15) is normalized by
%lambda_max*vmax. 19 - 25 normalized by K

lambdaMax = 1/(1/(sigma*n0));
maxTe = 25;
vMax = 1/(sqrt(q*maxTe/massIon));
DMax = lambdaMax*vMax;
kNorm = 1/(10e-18);
nNorm = 1/n0;
%Functions
F=[
(2*(x(22)/kNorm) + (x(21)/kNorm))*(x(4)/nNorm) * (x(1)/nNorm) - (k9+1/tau)*(x(2)/nNorm); %1, (2*k4+k3)*ne*nCl2-(k9+1/tau)nCl
(x(1)/nNorm)-n0*(1-delta) + 0.5*(x(2)/nNorm); %2, nCl2-n0(1-delta) + 0.5*nCl 
(x(3)/nNorm)- n0*delta;%3,%nAr- n0*delta
((x(20)/kNorm)+(x(21)/kNorm))*(x(1)/nNorm)*(x(4)/nNorm)-((x(7)/nNorm)+(x(8)/nNorm)+(x(6)/nNorm))*K*(x(5)/nNorm)-(x(25)/kNorm)*(x(5)/nNorm)*(x(4)/nNorm); %4, (k3+k3)*nCl2*ne-K(nCl_pos+nAr_pos+nCl2_pos)*nCl_neg-k8*nCl_neg*ne 
(x(19)/kNorm)*(x(1)/nNorm)*(x(4)/nNorm)-((x(14)/vMax)/x(13) + K*(x(5)/nNorm))*(x(6)/nNorm); % 5, k1*nCl2*ne-(v/dc + nCl-)*nCl2+
((x(23)/kNorm)*(x(2)/nNorm)+(x(20)/kNorm)*(x(1)/nNorm))*(x(4)/nNorm)-((x(14)/vMax)/x(13) + K*(x(5)/nNorm))*((x(7)/nNorm));%6 (k5*nCl+k2*nCl2)*ne-(v/dc+nCl_neg)*nCl_pos
(x(24)/kNorm)*(x(3)/nNorm)*(x(4)/nNorm)-((x(14)/vMax)/x(13)+ K*(x(5)/nNorm))*(x(8)/nNorm); %7 x(7, k6*nAr*ne-(v/dc  + k7*nCl_neg)*nAr_pos 
(x(5)/nNorm)+(x(4)/nNorm)-(x(6)/nNorm)-(x(7)/nNorm)-(x(8)/nNorm); %8 nCl_neg+ne-nCl2_pos+nCl_pos+nAr_pos
(x(14)/vMax) - sqrt(q/massIon*x(9))*sqrt(1+(x(11))/(1+(x(11)*x(12))));%9, v- sqrt((1+BetaS)/(1+BetaS*gamma_T))*sqrt(T)
x(10)-(x(5)/nNorm)/(x(4)/nNorm); %10 beta-nCl_neg/ne
x(10) - x(11)*exp((1+x(11))*(x(12)-1)/(2*(1+x(11)*x(12))));%11, Beta - BetaS*exp((1+BetaS)*(gamma_T-1)/(2*(1+BetaS*gamma_T)))
2*((x(19)/kNorm)*(x(1)/nNorm)+(x(23)/kNorm)*(x(2)/nNorm)+(x(24)/kNorm)*(x(3)/nNorm))+((x(20)/kNorm)-(x(21)/kNorm))*(x(1)/nNorm)-(1+x(10))*(2*((x(14)/vMax)/x(13)+K*x(10)*(x(4)/nNorm)))+(x(25)/kNorm)*(x(4)/nNorm)*x(10);%12;%17, 2*(k1*nCl2+k5*nCl+k6*nAr)+(k3-k3)*nCl2-(1+Beta)*(2*v/dc+k7*Beta*ne)+k8*ne*Beta
x(13)-0.5*R*L/(R*x(17)+L*x(18)); %13 dc-0.5*R*L/(R*hl+L*hr);
W-q*(x(14)/vMax)*7*x(9)*((x(6)/nNorm)+(x(7)/nNorm)+(x(8)/nNorm));%14 W-v*7x(9)(nCl2_pos+nCl_pos+nAr_pos)
(x(15)/DMax)-(x(16)/lambdaMax)*sqrt(q*x(9)/(2*massIon))*(1+x(12)+2*x(12)*x(10))/(1+x(10)*x(12));
(x(16)/lambdaMax)-1/(sigma*n0);
x(17)-(1+2*x(11)/x(12))/(1+x(11))*0.86/sqrt(3+L/((x(16)/lambdaMax)*2)+(.86*L*(x(14)/vMax)/(pi*x(12)*(x(15)/DMax)))^2);
x(18)-(1+3*x(11)/x(12))/(1+x(11))*0.8/sqrt(4+R/(x(16)/lambdaMax)+(.8*R*(x(14)/vMax)/(2.405*0.43*x(12)*(x(15)/DMax)))^2);
(x(19)/kNorm) - Act(1)/kNorm*exp(-B(1)/x(9));%k1
(x(20)/kNorm) - Act(2)/kNorm*exp(-B(2)/x(9));%k2
(x(21)/kNorm) - Act(3)/kNorm*exp(-B(3)/x(9));%k3
(x(22)/kNorm) - Act(4)/kNorm*exp(-B(4)/x(9));%k4
(x(23)/kNorm) - Act(5)/kNorm*exp(-B(5)/x(9));%k5
(x(24)/kNorm) - Act(6)/kNorm*exp(-B(6)/x(9));%k6
(x(25)/kNorm) - Act(7)/kNorm*exp(-B(7)/x(9));
x(9)/x(12)-TeV - (.5-TeV)/P];%k8


% 

end