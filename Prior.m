function F = Prior(current,i)
global center
global sd
global priorRecord
priorTime = tic;
prior = 0; %for when you are not doing cwise MH
prior = log(lognpdf(current(i),center(i),sd(i))) + prior;
F = prior;
priorRecord = [priorRecord prior];
priorTimeElapsed = toc(priorTime);
assignin('base', 'priorTimeElapsed', priorTimeElapsed);
end