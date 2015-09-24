function F = CalcEtchRate(i,x4)
global q_e
global ExpParameters
global TeSet;
nCl2_pos = x4(6);
nCl_pos = x4(7);
nAr_pos = x4(8);
nCl = x4(2);
Te = x4(9);
v = x4(12);
SCl = .3;
KCl = .5;
density_MgO = 3.56;
M = 40.3;
Na = 6.02e+23;
RxnProb = .25; %SCl*KCl;
FluxCl = nCl*v;
FluxAr_pos = nAr_pos*v;
TotalPosFlux = v*(nAr_pos+nCl_pos+ nCl2_pos);
Vdc = ExpParameters(i,6)*q_e;%-BiasFactor*Prf/TotalPosFlux;

%Surface Kinetics Problem
E = 6.0*Te + -Vdc; 
if E<0
    E = 0;
end

eos = 40; %eV
eod = 10; % eV
A = .05;
B = A;
Yd = B*(sqrt(E)-sqrt(eod));
Ys = A*(sqrt(E)-sqrt(eos));
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
EtchRate = EtchRate*6e+8*M/(Na*density_MgO); %nm/min
assignin('base', 'EtchRate', EtchRate);
F = EtchRate; 
end