function [ c ] = dotProduct( a, b, dim )
% DOT  Vector dot product.
%   C = DOT(A,B) returns the scalar product of the vectors A and B.
%   A and B must be vectors of the same length.  When A and B are both
%   column vectors, DOT(A,B) is the same as A'*B.
%
%   DOT(A,B), for N-D arrays A and B, returns the scalar product
%   along the first non-singleton dimension of A and B. A and B must
%   have the same size.
%
%   DOT(A,B,DIM) returns the scalar product of A and B in the
%   dimension DIM.

    if nargin == 2
        c = sum(a.*b);
    else
        c = sum(a.*b, dim);
    end
end