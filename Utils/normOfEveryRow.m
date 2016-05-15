function [ normVect ] = normOfEveryRow( matrix )
% normOfEveryRow: returns the 2-norm or Euclidean norm of vector v 
% where v is a row of a input matrix 
%%
    normVect = sqrt(sum(abs(matrix).^2, 2));
end