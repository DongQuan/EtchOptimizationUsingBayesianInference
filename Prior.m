function F = Prior(current,i)
global mean
global sd
priorTime = tic;
prior = 0;
prior = log(lognpdf(current(i),mean(i),sd(i))) + prior;
F = prior;
priorTimeElapsed = toc(priorTime);
assignin('base', 'priorTimeElapsed', priorTimeElapsed);
end