function F = plotter(theta)
global noUnknowns
global priorCenter
global priorSD
global proposalUB
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
end