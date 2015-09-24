function F = ProposeParameters(i)
global NoUnknowns;
global mean;
global sd;
Parameter = normrnd(mean(i),sd(i));
Parameter = abs(Parameter);
F = Parameter;
end

% %Non-dimensionalizations
% k7_nond = 5e-8;
% 
% %Test Values
% 
% mean(1) = 6e-10/k7_nond;
% mean(2) = 2.1e-10/k7_nond;
% mean(3) = 2e-12/k7_nond;
% mean(4) = 4e-9/k7_nond;
% mean(5) = 2e-10/k7_nond;
% mean(6) = 4e-11/k7_nond;
% mean(7) = 2.1e-10/k7_nond; % not exact values
% mean(8) = 10;
% for i=1:NoUnknowns
%     sd(i) = 100;
% end
