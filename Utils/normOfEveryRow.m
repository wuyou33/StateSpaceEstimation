function [ normVect ] = normOfEveryRow( matrix, direction )
    % normOfEveryRow. Returns the 2-norm or Euclidean norm of vector v, where v is a row of a input matrix.
    %
    %   [ normVect ] = normOfEveryRow( matrix, direction )
    %
    %   INPUT
    %       matrix      array of vectors;
    %       direction   direction.
    %
    %   OUTPUT
    %       normVect    array of Euclidean norm for every input vector.
    %
    narginchk(1, 2);
    
    if nargin == 1
        direction = 1;
    end
    
    normVect = sqrt(sum(abs(matrix).^2, direction));
end
