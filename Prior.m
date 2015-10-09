function F = Prior(current,i)
global mean
global sd
prior = 0;
prior = log(normpdf(current(i),mean(i),sd(i))) + prior;
F = prior;
end