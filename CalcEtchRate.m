function F = CalcEtchRate(plasmaVariables,expNo)
global expParameters %ion mass
global q
expNo = 31;

BiasFactor = expParameters(expNo,6);
nCl2_pos = plasmaVariables(6);
nCl_pos = plasmaVariables(7);
nAr_pos = plasmaVariables(8);
nCl = plasmaVariables(2);
Te = plasmaVariables(9);
v = plasmaVariables(14);
SCl = .3;
KCl = .5;
density_MgO = 3.56e+3; %kg/m^3
M = 40.3;
Na = 6.02e+23;
RxnProb = .25; %SCl*KCl;
FluxCl = nCl*v;
FluxAr_pos = nAr_pos*v;
TotalPosFlux = v*(nAr_pos+nCl_pos+ nCl2_pos);
Vdc = BiasFactor*q;%-BiasFactor*Prf/TotalPosFlux;

%Surface Kinetics Problem
E = 6.0*Te + -Vdc; 
if E<0
    E = 0;
end

eos = 40; %eV
eod = 10; % eV
A1 = .05;
B1 = A1;
Yd = B1*(sqrt(E)-sqrt(eod));
Ys = A1*(sqrt(E)-sqrt(eos));
if Ys<0
    Ys = 0;
end
if Yd<0
    Yd = 0;
end
assignin('base', 'E', E);
assignin('base', 'Ys', Ys);
assignin('base', 'Yd', Yd);
assignin('base', 'FluxAr_pos', FluxAr_pos);
assignin('base', 'FluxCl', FluxCl);
assignin('base', 'v', v);
assignin('base', 'TotalPosFlux', TotalPosFlux);
assignin('base', 'nCl_pos', nCl_pos);
assignin('base', 'FluxAr_pos', FluxAr_pos);
assignin('base', 'nAr_pos', nAr_pos);
assignin('base', 'nCl', nCl);
assignin('base', 'v', v);
assignin('base', 'Vdc', Vdc);
assignin('base', 'Te', Te);
EtchRate = (1-RxnProb*FluxCl/(RxnProb*FluxCl + Yd*TotalPosFlux))*(RxnProb*FluxCl+Ys*FluxAr_pos);
EtchRate = EtchRate*6e+7*M/(Na*density_MgO); %nm/min %6 comes from s to min
assignin('base', 'EtchRate', EtchRate);
F = EtchRate; 
end