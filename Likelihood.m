function [F] = Likelihood(current)
global data;
global likelihoodRecord
Like = 0;
likelihoodTime = tic;
noise = .1;
kGuess = testArr(current);
Like = normpdf(data,kGuess,noise) + 10e-20;
Like = log(Like) + Like;
% for i=1:length(data)
%     [etchGuess] = GlobalSolver(current,i);
%     Like = normpdf(data(i),etchGuess,current(15));
%     if (Like~=0)
%         Like = log(Like)+Like;
%     end
% end
likelihoodRecord = [likelihoodRecord Like];
likelihoodTimeElapsed = toc(likelihoodTime);
assignin('base', 'likelihoodTimeElapsed', likelihoodTimeElapsed);
F = Like;
