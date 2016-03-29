function [orderOfIndices] = prioritize(current)
stationary = current;
global center
global sd
global noUnknowns
global data
stationary = current;
uncertainty = zeros(noUnknowns,1);
for i = 1:noUnknowns
    for (j = 1:length(data))
        min = 10e-10; %put small min to not screw up noise
        max = center(i)+ 2*sd(i);
        stationary(i) = min;
        minK = testArr(current,j);
        stationary(i) = max;
        maxK = testArr(current,j);
        uncertainty(i) = abs(minK -maxK)+ uncertainty(i);
        stationary = current;
    end
end

[~,~,rnk] = unique(uncertainty);
orderOfIndices = rnk;

%how to rank in descending order?
%how to efficiently find index for j