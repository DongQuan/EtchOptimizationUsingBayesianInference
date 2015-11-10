function [F] = Likelihood(current)
global data;
global likelihoodRecord
Like = 0;
likelihoodTime = tic;
noise = 10; %change this back to unknown noise parameter eventually

for i=1:length(data)
    [etchGuess] = GlobalSolver(current,i);
    Like = normpdf(data(i),etchGuess,noise) + 10e-20;
    Like = log(Like) + Like;
end
%likelihoodRecord = [likelihoodRecord Like];
likelihoodTimeElapsed = toc(likelihoodTime);
assignin('base', 'likelihoodTimeElapsed', likelihoodTimeElapsed);
F = Like;
