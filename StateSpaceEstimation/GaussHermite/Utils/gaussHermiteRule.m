function [ points, weights ] = gaussHermiteRule( order, dim )
    % gaussHermiteQuadratureRule. Draw cubature points and corresponding weights using Gauss-Hermite cubature rules.
    %
    %   Generate points and weights according to Gauss-Hermite quadrature rule.
    %   Weight function is chosen to be the standard Gaussian density with zero mean and unit variance N(0; I).
    %   The interval of interest is chosen to be (-infinity; +infinity).
    %   According to the fundamental theorem of Gauss-Hermite quadrature, the quadrature points are chosen to be the zeros of the m-th order Hermite polynomial.
    %   Since the zeros of the Hermite polynomials are distinct,
    %   it is noteworthy that the determinant of the coefficient matrix in is the well known Vandermonde's determinant that is nonzero.
    %   For an m-point quadrature scheme, the resulting quadrature rule is exact for all polynomials of degree  <= 2m - 1
    %
    %   Suppose J is a symmetric tridiagonal matrix with zero diagonal elements and J(i, i+1) = sqrt(i/2), 1 <= i <= m-1
    %   Then the quadrature point KSI is taken to be KSI(k) = sqrt(2)*eigenvalue_J(k);
    %   where "l is the l-th eigenvalue of J;
    %   and the corresponding weight w(k) = v(k, 1)^2 where v(k, 1) is the first element of the k-th normalized eigenvector of J.
    %
    %   [ points, weights ] = gaussHermiteRule( order, dim )
    %
    %   INPUT
    %       order   order of Gauss-Hermite polynomial;
    %       dim   	dimension.
    %
    %   OUPTUT
    %       points  	generated points;
    %       weights 	generated weights for every point.
    %
    narginchk(1, 2);
    
    indices   = 1 : order - 1;
    elements  = sqrt(indices / 2);
    stdm = diag(elements, 1) + diag(elements, -1);
    
    [ eigVal, eigVect ] = eig(stdm);
    [points, ind] = sort(diag(eigVect));
    
    eigVal   = eigVal(:, ind)';
    weights  = sqrt(pi) * eigVal(:, 1).^2;
    
    points = points'*sqrt(2);
    weights = weights'/sum(weights);
    
    if nargin == 2
        replRoots(1:dim) = {points};
        points = cartesianProduct(replRoots)';
        
        replWeights(1:dim) = {weights};
        tmpWeights = cartesianProduct(replWeights)';
        weightLenght = size(tmpWeights, 2);
        weights(1:weightLenght) = prod(tmpWeights(:, 1:weightLenght));
    end
end
