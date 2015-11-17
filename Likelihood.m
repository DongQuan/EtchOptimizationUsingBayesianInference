function [F] = Likelihood(current)
global data;
global likelihoodRecord
global etchRecord
likeData = 0;
likelihoodTime = tic;
noise = current(15); %change this back to unknown noise parameter eventually

for i=1:length(data)
    [etchGuess] = GlobalSolver(current,i);
    etchRecord = [etchRecord etchGuess]
    Like = normpdf(data(i),etchGuess,noise) + 10e-20;
    likeData = log(Like) + likeData;
end
likelihoodRecord = [likelihoodRecord likeData];
likelihoodTimeElapsed = toc(likelihoodTime);
assignin('base', 'likelihoodTimeElapsed', likelihoodTimeElapsed);
F = likeData;
