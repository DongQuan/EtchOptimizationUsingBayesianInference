function F = FactorialDesign()
% Generate factorial design
Pressure = [.5 2]
flowRate = [10 40]
ICPPower = [300 700]
RFPower = [100 200]
delta = [0 .5 1]
dcbias = [0 300]
ind =1;
Exp = zeros(length(Pressure)*length(flowRate)*length(ICPPower)*length(RFPower)*length(delta)*length(dcbias),6);
for P=1:length(Pressure)
    for f=1:length(flowRate)
        for ICP = 1:length(ICPPower)
            for RF = 1:length(RFPower)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
                 for d = 1:length(delta)
                       for dc= 1:length(dcbias)
                           Exp(ind,5) = delta(d)     
                           Exp(ind,4) = RFPower(RF)       
                           Exp(ind,3) = ICPPower(ICP) 
                           Exp(ind,2) = flowRate(f);
                            Exp(ind,1) = Pressure(P)
                            Exp(ind,6) = dcbias(dc)  
                            ind = ind+1;
                       end
                 end
            end
        end
    end
end
F = Exp
end