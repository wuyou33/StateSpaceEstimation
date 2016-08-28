function [ c ] = dotProduct( a, b, dim )
    % DOTPRODUCT Vector dot product.
    %   c = dotProduct(a, b) returns the scalar product of the vectors A and B.
    %   a and b must be vectors of the same length. When a and b are both.
    %
    %   WORKS ONLY WITH REAL NUMBERS.
    %
    %   dotProduct(a, b), for N-D arrays a and b, returns the scalar product
    %   along the first non-singleton dimension of a and b. a and b must have the same size.
    %
    %   dotProduct(a, b, dim) returns the scalar product of a and b in the dimension dim.
    
    if nargin == 2
        c = sum(a.*b);
    else
        c = sum(a.*b, dim);
    end
end
