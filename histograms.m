clear
close all
global noUnknowns
%Make histograms of parameter estimates

% x = (0:.00001:10)';
% 
% for (i=1:100)
%     test(i) = lognrnd(log(2.1e-18),10e-10);
% end
% y = lognpdf(x,log(2.1e-18),10e-18);
% %y2 = lognpdf(x,(2.1e-18),1e-18);
% figure(1);
% plot(x,y)
% %figure(2) 
% %plot(x,y2)

theta = xlsread('SyntheticData.xlsx','theta-50','A1:O750');
for (i=1:noUnknowns)
    subplot(3,5,i)
    hist(theta(:,i)); 
    xlabel(i)
end
kmeans = mean(theta);
k7_nond = 5e-8;
k1 = 3e-16;
k2 = 2.1e-18;
k3 = 1.5e-16;
k4 = 3e-15;
k5 = 1e-16;
k6 = 2e-17;
k8 = k2;
RealValues = [3e-16 2.1e-18 1.5e-16 3e-15 1e-16 2e-17 k2 2 5 6 10 18 5 8 0];
Diff = (kmeans - RealValues)./RealValues;