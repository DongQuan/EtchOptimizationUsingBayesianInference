function F = testMatrix(current,i)
global epMatrix
result = current*epMatrix(:,i);
F = result;
end
