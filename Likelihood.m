function [F] = Likelihood(current)
global data;
Like = 0;
likelihoodTime = tic;
for i=1:length(data)
    [etchGuess] = GlobalSolver(current,i)
    Like = normpdf(data(i),etchGuess,current(15));
    if (Like~=0)
        Like = log(Like)+Like;
    end
end
likelihoodTimeElapsed = toc(likelihoodTime);
assignin('base', 'likelihoodTimeElapsed', likelihoodTimeElapsed);
F = Like
