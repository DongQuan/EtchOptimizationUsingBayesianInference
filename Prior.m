function F = Prior(current,subBlock)
global center
global sd
global priorRecord
priorTime = tic;
prior = 0; %for when you are not doing cwise MH
kNorm = 10e+18;
if subBlock == 1
    for index = 1:7
        prior = log(lognpdf(current(index)*kNorm,center(index),sd(index))) + prior;
    end
elseif subBlock == 2
    for actIndex = 8:14
        prior = log(lognpdf(current(actIndex),center(actIndex),sd(actIndex))) + prior;
    end
else
        prior = log(lognpdf(current(15),center(15),sd(15)));
end
F = prior;
priorRecord = [priorRecord prior];
priorTimeElapsed = toc(priorTime);
assignin('base', 'priorTimeElapsed', priorTimeElapsed);
end