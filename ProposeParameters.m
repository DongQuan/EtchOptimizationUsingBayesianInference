function F = ProposeParameters(i)
global center
global sd
global proposedParameterRecord
parameter = lognrnd(center(i),sd(i));
%parameter = lognrnd(mean(i),sd(i));
% kNorm = 10e+18;
% %scale prexponential factor values back down
% if i<=7
%     parameter = parameter/kNorm;
% end
% %parameter = abs(parameter);
proposedParameterRecord = [proposedParameterRecord parameter];
F = parameter;
end
