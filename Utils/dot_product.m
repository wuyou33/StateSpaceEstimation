function [ c ] = dot_product(a, b, dim)
    % dot_product. Vector dot product.
    %
    %   [ c ] = dot_product( a, b, dim )
    %
    %   Returns the scalar product of the vectors a and b. a and b must be vectors of the same length. When a and b are both.
    %   WORKS ONLY WITH REAL NUMBERS.
    %
    %   for N-D arrays a and b, returns the scalar productalong the first non-singleton dimension of a and b. a and b must have the same size.
    %
    %   returns the scalar product of a and b in the dimension dim.
    %
    %   INPUT
    %       a       first vector;
    %       b       second vector;
    %       dim     dimension.
    %
    %   OUTPUT
    %       c   result of scalar product.
    %
    if nargin == 2
        c = sum(a.*b);
    else
        c = sum(a.*b, dim);
    end
end
