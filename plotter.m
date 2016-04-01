function F = plotter(theta,N,trainingExp,testExp,trainingData,real)
global noUnknowns
global priorCenter
global priorSD
global proposalUB
global proposalLB
global expParameters
global data
draws = N/2;
figure;
    for i =1:noUnknowns
        subplot(3,3,i);
        outputTitle = sprintf('Unknown %d',i);
        paramEdges = [0:2:proposalUB(i)];
        f = histogram(theta(:,i),paramEdges);
        f.Normalization = 'probability';
        counts = f.Values;
        hold on
        line([real(i) real(i)],[0 max(counts)], 'Color', 'r')
        hold on
        xmin = min(theta(:,i));
        xmax = max(theta(:,i));
        x = xmin:1:xmax;
        lnorm = lognpdf(x,priorCenter(i),priorSD(i));
        plot(x,lnorm,'g');
        title(outputTitle);
        xlabel('Value');
        ylabel('Frequency');
    end
param = zeros(1,noUnknowns);
predTraining = zeros(draws,length(data));
predTest = zeros(draws,length(testExp));


for j = 1:draws
    for p = 1:noUnknowns
        param(p) = theta(round(rand(1)*(N-1))+1,p);
    end
    for i =1:length(data)
        predTraining(j,i) = testArr(param,i);
    end
end
predTraining

figure; %Training Experiments
for k =1:length(trainingExp)
    outputTitle = sprintf('Train Exp# %d',k);
    edges = [round((min(predTraining(:,k))-100)):10:round((max(predTraining(:,k))+100))]
    subplot(2,3,k);
    histogram(predTraining(:,k),edges);
    [p,edges] = histcounts(predTraining(:,k),edges);
    hold on
    line([data(k) data(k)],[0 max(p)], 'Color', 'r')
    title(outputTitle);
    xlabel('Value');
    ylabel('Frequency');
end

expParameters = testExp;
for j = 1:draws
    for p=1:noUnknowns
        param(p) = theta(round(rand(1)*(N-1))+1,p);
    end
    for i =1:length(testExp)
        predTest(j,i) = testArr(param,i);
    end
end

for i=1:length(testExp)
    testData(i) = testArr(real,i); 
end

figure; %Test Experiments

for k =1:length(testExp)
    outputTitle = sprintf(' Test Exp# %d',k);
    edges = [round((min(predTest(:,k))-100)):10:round((max(predTest(:,k))+100))];
    subplot(3,3,k);
    histogram(predTest(:,k),edges);
    [p,edges] = histcounts(predTest(:,k),edges);
    hold on
    line([testData(k) testData(k)],[0 max(p)], 'Color', 'r')
    title(outputTitle);
    xlabel('Value');
    ylabel('Frequency');
end

figure; 
for k =1:noUnknowns
    outputTitle = sprintf('Unknown %d',k);
    subplot(3,3,k);
    plot(theta(:,k));
    hold on
    line([0 N],[real(k) real(k)], 'Color', 'r')
    hold on
    line([0 N],[proposalLB(i) proposalLB(i)], 'Color', 'g')
    line([0 N],[proposalUB(i) proposalUB(i)], 'Color', 'g')
    title(outputTitle);
    xlabel('Cycle #');
    ylabel('Value');
end
end