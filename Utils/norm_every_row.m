function [ norm_vect ] = norm_every_row(matrix, direction)
    % norm_every_row. Returns the 2-norm or Euclidean norm of vector v, where v is a row of a input matrix.
    %
    %   [ norm_vect ] = norm_every_row( matrix, direction )
    %
    %   INPUT
    %       matrix      array of vectors;
    %       direction   direction.
    %
    %   OUTPUT
    %       norm_vect    array of Euclidean norm for every input vector.
    %
    narginchk(1, 2);
    
    if nargin == 1
        direction = 1;
    end
    
    norm_vect = sqrt(sum(abs(matrix).^2, direction));
end
