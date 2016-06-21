function [F] = Likelihood(current)
global data;
global likelihoodRecord
global etchRecord
imageRows = 10;
imageColumns = 10;
likeData = 0;
noise = 0.1;%current(end); %change this back to unknown noise parameter eventually
for i=1:length(data)
    [etchGuess,xmulti,flag] = GlobalSolver(current,i);
    etchRecord = [etchRecord etchGuess];
    Like = normpdf(data(i),etchGuess,noise) + 10e-20;
    likeData = log(Like) + likeData;
end
likelihoodRecord = [likelihoodRecord likeData];
F = likeData;
