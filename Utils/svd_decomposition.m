function [ sqrt_matrix ] = svd_decomposition(matrix)
    % svd_decomposition. Calculate square root factor of matrix through Singular value decomposition.
    %
    %   Produce square root decompoistion from standard svd decomposition via following equation
    %
    %       sqrtDecomposition = 0.5*(u + v) * sqrt(s).
    %
    %   Another similar equation:
    %       sqrtDecomposition = u * sqrt(s)
    %
    %   Where u, s, v defined as:
    %   [u, s, v] = svd(x) produces a diagonal matrix s, of the same dimension as x and with nonnegative diagonal elements in
    %   decreasing order, and unitary matrices U and V so, that X = U*S*V'.
    %
    %   More details are here:
    %
    %       http://stats.stackexchange.com/questions/238963/how-to-do-svd-instead-of-cholesky-for-ltl
    %
    %   INPUT
    %       matrix  input matrix.
    %
    %   OUPUT
    %       sqrt_matrix    square root decomposition of input matrix (0.5*(u + v) * sqrt(s), where u, v, s - results of svd decomponsition).
    %
    [u, s, v] =  svd(matrix);
    sqrt_matrix = 0.5*(u + v) * sqrt(s);
end
