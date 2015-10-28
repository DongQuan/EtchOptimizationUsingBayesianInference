function F = FactorialDesign()
global q
global massIon
% Generate factorial design
Pressure = [.5 1.5];
flowRate = [10 40];
ICPPower = [400 700];
RFPower = [100 200];
delta = [0 .5 1];
dcbias = [0];
ind =1;
ExpSet = zeros(length(Pressure)*length(flowRate)*length(ICPPower)*length(RFPower)*length(delta)*length(dcbias),9);
for P=1:length(Pressure)
    for f=1:length(flowRate)
        for ICP = 1:length(ICPPower)
            for RF = 1:length(RFPower)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
                 for d = 1:length(delta)
                       for dc= 1:length(dcbias)
                           Te = delta(d)/(Pressure(P))+(ICPPower(ICP)+RFPower(RF))/1000 +2;
                           nCl = -(delta(d)-1)^2+1.2;
                           v = sqrt(q*Te/massIon);
                           ExpSet(ind,5) = delta(d);    
                           ExpSet(ind,4) = RFPower(RF);       
                           ExpSet(ind,3) = ICPPower(ICP) ;
                            ExpSet(ind,2) = flowRate(f);
                            ExpSet(ind,1) = Pressure(P);
                            ExpSet(ind,6) = dcbias(dc) ; 
                            ExpSet(ind,7) = Te;
                            ExpSet(ind,8) = nCl;
                            ExpSet(ind,9) = v;
                            ind = ind+1;
                       end
                 end
            end
        end
    end
end
F = ExpSet;
end